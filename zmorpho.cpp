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
#include <algorithm>

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

/// Main procedure for curvature microscope.
int main(int argc, char** argv) {
    unsigned int zoom=2;
    std::string inLL, outLL, sOutImage;
    CmdLine cmd; cmd.prefixDoc = "\t";
    cmd.add( make_option('z',zoom).doc("Integer zoom factor") );

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

    size_t w, h;
    unsigned char* inIm = io_png_read_u8_gray(argv[1], &w, &h);
    if(! inIm) {
        std::cerr << "Error reading as PNG image: " << argv[1] << std::endl;
        return 1;
    }

    Timer timer;

    Rect rectSelect= {0,0,(int)w,(int)h};
    const Rect R={int(zoom*MARGIN), int(zoom*MARGIN),
                  int(zoom*rectSelect.w), int(zoom*rectSelect.h)}; // Image ROI
    rectSelect.x -= MARGIN;
    rectSelect.y -= MARGIN;
    rectSelect.w += 2*MARGIN;
    rectSelect.h += 2*MARGIN;

    unsigned char* inImage = extract(inIm,w,h, rectSelect);
    free(inIm);

    std::cout << "Tree extraction. " << std::flush;
    timer.tick();
    LsTree tree(inImage, rectSelect.w, rectSelect.h);
    std::cout << "Shapes: " << tree.iNbShapes << ". " << std::flush;
    timer.time();
    delete [] inImage;

    std::cout << "Level lines smoothing. " << std::flush;
    timer.tick();
    float scale = zoom*.5f;
    LLTree* lltree = convert(tree,zoom);
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
        for(std::vector<DPoint>::iterator it=dline.begin(); it!=dline.end(); ++it)
            line.push_back( Point((float)it->x,(float)it->y) );
    }
    timer.time();

    std::cout << "Image reconstruction. " << std::flush;
    timer.tick();
    inImage = reconstruct(*lltree, zoom*rectSelect.w, zoom*rectSelect.h, R);
    if(io_png_write_u8(argv[2], inImage, R.w, R.h, 1)!=0) {
        std::cerr << "Error writing image file " << sOutImage << std::endl;
        return 1;
    }
    delete [] inImage;
    timer.time();

    delete [] lltree->root()->ll;
    delete lltree;

    return 0;
}
