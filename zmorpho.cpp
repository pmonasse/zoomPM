// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file zmorpho.cpp
 * @brief Morphological zoom
 * 
 * Copyright (c) 2023 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "boundaries/tree.h"
#include "microCurv/lltree.h"
#include "microCurv/gass.h"
#include "microCurv/fill_curve.h"
#include "microCurv/image.h"
#include "io_png.h"
#include "cmdLine.h"
#include "xmtime.h"
#include <cstring>

/// Timer class to measure real time (not CPU time)
class Timer {
    unsigned long t; ///< Current time in milliseconds
public:
    Timer() { tick(); } ///< Constructor
    unsigned long tick() { return t=xmtime(); } ///< Reset time
    void time() { ///< Display elapsed time and reset current time
        unsigned long told = t;
        std::cout << "Time = " << (tick()-told)/1000.0f << "s" << std::endl;
    }
};

/// Test whether the 3-channel image \a data is actually gray-scale.
bool all_gray(const unsigned char* data, int w, int h) {
    int n=w*h;
    for(int i=1; i<3; i++) {
        const unsigned char *p=data, *q=p+i*n;
        for(int j=0; j<n; j++)
            if(*p++ != *q++)
                return false;
    }
    return true;
}

/// The image is enlarged by this number of pixels on each side before
/// extraction of level lines, in order to reduce border effects.
static const int MARGIN=1;

LLTree* convert(const LsTree& t, unsigned int z) {
    const int n=t.iNbShapes;
    LevelLine* lls = new LevelLine[n];
    for(int i=0; i<n; i++) {
        const std::vector<LsPoint>& c = t.shapes[i].contour;
        LevelLine& ll=lls[i];
        ll.level = (float)t.shapes[i].gray;
        for(size_t j=0; j<c.size(); j++)
            ll.line.emplace_back(z*c[j].x-.5f,z*c[j].y-.5f);
        ll.line.push_back(ll.line.front()); // Close level line
    }
    std::vector<LLTree::Node> nodes;
    for(int i=0; i<n; i++)
        nodes.emplace_back(&lls[i]);
    LLTree* lltree = new LLTree(nodes);
    std::vector<LLTree::Node>& nds = lltree->nodes();
    for(int i=0; i<n; i++) {
        const LsShape& s = t.shapes[i];
        nds[i].parent = s.parent? &nds[0]+(s.parent-t.shapes):0;
        nds[i].sibling = s.sibling? &nds[0]+(s.sibling-t.shapes):0;
        nds[i].child = s.child? &nds[0]+(s.child-t.shapes):0;
    }
    return lltree;
}

/// Reconstruct image from level sets stored in \a tree.
unsigned char* reconstruct(LLTree& tree, int w, int h, const Rect& R) {
    unsigned char* outImage = new unsigned char[w*h];
    std::fill_n(outImage, w*h, 0);
    std::vector< std::vector<float> > inter;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        fill_curve(it->ll->line,(unsigned char)it->ll->level,
                   outImage,w,h, &inter);
    unsigned char* out = crop(outImage,w,h, R);
    delete [] outImage;
    return out;
}

unsigned char* zoom_channel(const unsigned char* inIm, int w, int h, int z) {
    Timer T;

    Rect rectSelect= {0,0,(int)w,(int)h};
    const Rect R={int(z*MARGIN), int(z*MARGIN),
                  int(z*rectSelect.w), int(z*rectSelect.h)}; // Image ROI
    rectSelect.x -= MARGIN;
    rectSelect.y -= MARGIN;
    rectSelect.w += 2*MARGIN;
    rectSelect.h += 2*MARGIN;

    unsigned char* inImage = extract(inIm,w,h, rectSelect);

    std::cout << "Tree extraction. " << std::flush;
    T.tick();
    LsTree tree(inImage, rectSelect.w, rectSelect.h);
    std::cout << "Shapes: " << tree.iNbShapes << ". " << std::flush;
    T.time();
    delete [] inImage;

    std::cout << "Level lines smoothing. " << std::flush;
    T.tick();
    float scale = z*.5f;
    LLTree* lltree = convert(tree,z);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<tree.iNbShapes; i++) {
        std::vector<Point>& line = lltree->nodes()[i].ll->line;
        std::vector<DPoint> dline;
        for(std::vector<Point>::iterator it=line.begin(); it!=line.end(); ++it)
            dline.push_back( DPoint((double)it->x,(double)it->y) );
        assert(dline.front()==dline.back());
        gass(dline, 0.0, scale);
        line.clear();
        for(std::vector<DPoint>::iterator it=dline.begin();it!=dline.end();++it)
            line.push_back( Point((float)it->x,(float)it->y) );
    }
    T.time();

    std::cout << "Image reconstruction. " << std::flush;
    T.tick();
    inImage = reconstruct(*lltree, z*rectSelect.w, z*rectSelect.h, R);
    T.time();

    delete [] lltree->root()->ll;
    delete lltree;
    return inImage;
}

/// Main procedure for curvature microscope.
int main(int argc, char** argv) {
    unsigned int zoom=2;
    bool bColor=false;
    std::string inLL, outLL, sOutImage;
    CmdLine cmd; cmd.prefixDoc = "\t";
    cmd.add( make_option('z',zoom).doc("Integer zoom factor") );
    cmd.add( make_option('c',bColor,"color").doc("Handle color image") );

    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << "Error: " << s << std::endl;
        return 1;
    }
    if(argc!=3) {
        std::cerr << "Usage: " <<argv[0]<< " [options] imgIn.png imgOut.png\n"
                  << "Options:\n" << cmd << std::endl;
        return 1;
    }
    if(zoom==0) {
        std::cerr << "The value of z must be positive" << std::endl;
        return 1;
    }

    size_t w, h, channels=1;
    unsigned char* inIm;
    if(bColor) {
        inIm = io_png_read_u8_rgb(argv[1], &w, &h);
        if(all_gray(inIm,w,h))
            std::cout << "Image is grayscale: process one channel" << std::endl;
        else
            channels=3;
    } else
        inIm = io_png_read_u8_gray(argv[1], &w, &h);
    if(! inIm) {
        std::cerr << "Error reading as PNG image: " << argv[1] << std::endl;
        return 1;
    }

    unsigned char* out[3];
    for(int i=0; i<channels; i++) {
        if(channels>1)
            std::cout << "*Processing channel " << i << "*" << std::endl;
        out[i] = zoom_channel(inIm+i*w*h, w, h, zoom);
    }

    bool ok=true;
    if(channels==1)
        ok = (io_png_write_u8(argv[2], out[0], zoom*w, zoom*h, 1)==0);
    else {
        int sz = zoom*w*zoom*h;
        unsigned char* outImage = new unsigned char[sz*channels];
        for(int i=0; i<channels; i++)
            memcpy(outImage+sz*i, out[i], sz*sizeof(unsigned char));
        ok = (io_png_write_u8(argv[2], outImage, zoom*w, zoom*h, channels)==0);
        delete [] outImage;
    }
    for(int i=0; i<channels; i++)
        delete [] out[i];

    if(! ok) {
        std::cerr << "Error writing image file " << sOutImage << std::endl;
        return 1;
    }
    
    free(inIm);

    return 0;
}
