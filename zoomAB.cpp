// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file zoomAB.cpp
 * @brief Super-resolution with Perona-Malik diffusion (Belahmidi's algorithm)
 * 
 * Copyright (c) 2004 Abdelmounim Belahmidi
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

#include "cmdLine.h"
#include "io_png.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include "xmtime.h"

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
bool all_gray(const float* data, int w, int h) {
    int n=w*h;
    for(int i=1; i<3; i++) {
        const float *p=data, *q=p+i*n;
        for(int j=0; j<n; j++)
            if(*p++ != *q++)
                return false;
    }
    return true;
}

const float sqrt2=M_SQRT1_2;
const float sqrt22=sqrt2/2;
const float minGrad=0.001f; ///< Minimum gradient norm for applying anisotropy
const float dt = 0.1f; ///< Time step for diffusion

/// Compute the 0-order spline interpolation of the input image
void zoom_duplication(const float* in, int w, int h, int z, float* out) {
    const int zw=z*w;
    for (int j=0; j<h; j++) {
        const float* line=out;
        for(int i=0; i<w; i++, in++)
            for(int k=0; k<z; k++)
                *out++ = *in;
        for(int k=1; k<z; k++, out+=zw) // Duplicate z-1 times each line
            memcpy(out, line, zw*sizeof(float));
    }
}

/// Replace each block of \a z x \a z pixels in \a in by its average.
/// Write the result in \a average.
void projectionPU(const float*  in, int w, int h, int z,
                  float* temp, float* average) {
    float norm = 1.0f/float(z*z);
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int offset = w*j+i;
            float s=0.0;
            for(int y=j*z; y<(j+1)*z; y++)
                for(int x=i*z; x<(i+1)*z; x++){
                    int Zoffset = z*w*y+x;
                    s += in[Zoffset];
                }
            s *= norm;
            temp[offset]=s;
        }
    zoom_duplication(temp,w,h, z, average);
}

void derivatives(const float a0[3], const float a1[3], const float a2[3],
                 float& grad, float& g, float& h) {
    static const float c2=0.707106781186547129*4.0+4.0;
    float c= a2[2] - a0[0];
    float d= a2[0] - a0[2];
    float ax= a1[2] - a1[0] + sqrt22*(c-d);
    float ay= a2[1] - a0[1] + sqrt22*(c+d);
    float az=ax*ay;
    ax *= ax;
    ay *= ay;
    grad = ax+ay;
    if(grad < ::minGrad) {
        h = g = ((a0[0]+a0[2]+a2[0]+a2[2])*sqrt2 +
                 (a0[1]+a1[0]+a1[2]+a2[1]) - c2*a1[1])/c2;
        return;
    }
    float li=1.0/grad;
    az*=li;
    ax*=li;
    ay*=li;  
    li=az*az;
    float l0=-2.0+4.0*li;
    float l1=ay*(ay-ax);
    float l2=ax*(ax-ay);
    float l3=li-.5*az;
    float l4=l3+az;
    g = l0*a1[1] +
        l1*(a1[0]+a1[2]) + l2*(a0[1]+a2[1]) +
        l3*(a2[2]+a0[0]) + l4*(a0[2]+a2[0]); 
    h = l0*a1[1] +
        l2*(a1[0]+a1[2]) + l1*(a0[1]+a2[1]) +
        l4*(a2[2]+a0[0]) + l3*(a0[2]+a2[0]); 
}

void compute_derivatives(const float *im, int w, int dy,
                         float *grad, float *evolxixi, float *evolnn) {
    // Inside
    const int jmax = dy-1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j=1; j<jmax; j++) {
        const float *a0=im+(j-1)*w, *a1=a0+w, *a2=a1+w;     
        float *gg=grad+j*w+1, *g=evolxixi+j*w+1, *h=evolnn+j*w+1;
        for(int i=1; i+1<w; i++) {
            derivatives(a0,a1,a2, *gg, *g, *h);
            a0++; a1++; a2++; gg++; g++; h++;
        }
    }
    // Top boundary
    const float *a0=im, *a1=im+w;
    float *gg=grad, *g=evolxixi, *h=evolnn;
    for(int k=1; k+1<w; k++) {
        derivatives(a0, a0, a1, gg[k], g[k], h[k]);
        a0++; a1++;
    }
    // Bottom boundary
    a0=im+w*(dy-2); a1=a0+w;
    gg=grad+w*(dy-1); g=evolxixi+w*(dy-1); h=evolnn+w*(dy-1);
    for(int k=1; k+1<w; k++) {
        derivatives(a0, a1, a1, gg[k], g[k], h[k]);
        a0++; a1++;
    }
    // Left boundary
    gg=grad+w; g=evolxixi+w; h=evolnn+w;
    for(int k=1; k+1<dy; k++) {
        float t0[3]={im[0+(k-1)*w], im[0+(k-1)*w], im[1+(k-1)*w]};
        float t1[3]={im[0+(k+0)*w], im[0+(k+0)*w], im[1+(k+0)*w]};
        float t2[3]={im[0+(k+1)*w], im[0+(k+1)*w], im[1+(k+1)*w]};
        derivatives(t0, t1, t2, *gg, *g, *h);
        gg+=w; g+=w; h+=w;
    }
    // Right boundary
    gg=grad+2*w-1; g=evolxixi+2*w-1; h=evolnn+2*w-1;
    for(int k=1; k+1<dy; k++) {
        float t0[3]={im[w-2+(k-1)*w], im[w-1+(k-1)*w], im[w-1+(k-1)*w]};
        float t1[3]={im[w-2+(k+0)*w], im[w-1+(k+0)*w], im[w-1+(k+0)*w]};
        float t2[3]={im[w-2+(k+1)*w], im[w-1+(k+1)*w], im[w-1+(k+1)*w]};
        derivatives(t0, t1, t2, *gg, *g, *h);
        gg+=w; g+=w; h+=w;
    }
    // Corners
    { // Top-left
        float t0[3]={im[0], im[0], im[1]};
        float t1[3]={im[w],im[w],im[w+1]};
        derivatives(t0, t0, t1, grad[0], evolxixi[0], evolnn[0]);
    }
    { // Top-right
        int i = w-1;
        float t0[3]={im[i-1],   im[i],   im[i]};
        float t1[3]={im[i-1+w],im[i+w],im[i+w]};
        derivatives(t0, t0, t1, grad[i], evolxixi[i], evolnn[i]);
    }
    { // Bottom-left
        int i = (dy-1)*w;
        float t0[3]={im[i-w],im[i-w],im[i-w+1]};
        float t1[3]={im[i],   im[i],   im[i+1]};
        derivatives(t0, t1, t1, grad[i], evolxixi[i], evolnn[i]);
    }
    { // Bottom-right
        int i = w*dy-1;
        float t0[3]={im[i-w-1],im[i-w],im[i-w]};
        float t1[3]={im[i-1],   im[i],   im[i]};
        derivatives(t0, t1, t1, grad[i], evolxixi[i], evolnn[i]);
    }
}

/// Zoom of low-res image \a lr into high-res image \a hr.
void zoomAB(const float* lr, int w, int h, int z, float* hr) {
    const int zw=z*w, zh=z*h;
    const int iSizeImageZoom = zw*zh;
    float* temp = new float[w*h];
    float* dup     = new float[iSizeImageZoom];
    float* grad    = new float[iSizeImageZoom];
    float* uxixi   = new float[iSizeImageZoom];
    float* unn     = new float[iSizeImageZoom];
    float* average = new float[iSizeImageZoom];

    const int n_iter = int(z*z/::dt+0.5f);
    std::cout << "Iterations: " << n_iter << "." << std::endl;
    

    // Generate dup, the 0-order spline interpolation of the input image
    zoom_duplication(lr,w,h, z, dup);
    // Copy it as initial zoomed image hr
    memcpy(hr, dup, iSizeImageZoom*sizeof(float));

    for(int k=0; k<n_iter; k++) {
        compute_derivatives(hr,zw,zh, grad, uxixi, unn);
        projectionPU(hr,w,h, z, temp, average);

        const float* grd=grad;
        const float* uxi=uxixi;
        const float* un=unn;
        const float* avg=average;
        const float* uo=dup;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<iSizeImageZoom; i++) {
            float reac = average[i] - dup[i];
            float evolx = uxixi[i];
            float evoln = unn[i] / (1.0f + 0.01f * grad[i]);
            hr[i] += ::dt*(evolx + evoln - reac);
            if(hr[i] > 255.0)
                hr[i] = 255;
            else if(hr[i] <0.0)
                hr[i] = 0;
        }
    }
    delete [] temp;
    delete [] dup;
    delete [] uxixi;
    delete [] unn;
    delete [] grad;
    delete [] average;
}

int main(int argc, char* argv[]) {
    CmdLine cmd;
    bool bColor=false;
    unsigned int z=2; // Max number of pixels to discard
    cmd.add( make_option('z',z).doc("Integer zoom factor") );
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
    if(z==0) {
        std::cerr << "The value of z must be positive" << std::endl;
        return 1;
    }

    size_t w, h, channels=1;
    float* data;
    if(bColor) {
        data = io_png_read_f32_rgb(argv[1], &w, &h);
        if(all_gray(data,w,h))
            std::cout << "Image is grayscale: process one channel" << std::endl;
        else
            channels=3;
    } else
        data = io_png_read_f32_gray(argv[1], &w, &h);
    if(! data) {
        std::cerr << "Error loading image " << argv[1] << std::endl;
        return 1;
    }

    float* out = new float[channels*z*z*w*h];
    Timer T;
    for(int i=0; i<channels; i++) {
        if(channels>1)
            std::cout << "Channel " << i << ". " << std::flush;
        zoomAB(data+w*h*i, (int)w, (int)h, (int)z, out+z*z*w*h*i);
    }
    T.time();
    if(io_png_write_f32(argv[2], out, z*w, z*h, channels) != 0) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }

    free(data);
    delete [] out;
    return 0;
}