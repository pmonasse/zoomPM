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

const float minGrad=0.001f; ///< Minimum gradient norm for applying anisotropy
const float dt = 0.1f; ///< Time step for diffusion

/// Clamp \a f in interval [0,255].
void clamp(float& f) {
    if(f<0.0f) f=0.0f;
    else if(f>255.0f) f=255.0f;
}

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

/// Image \a in is of actual size (zw+2)x(zh+2) because it has a padding of
/// 1 pixel each side: fill this padding.
void fill_padding(float* in, int zw, int zh) {
    for(int i=1; i<=zh; i++) {
        int ii = i*(zw+2);
        in[ii]      = in[ii+1];  // Pad left
        in[ii+zw+1] = in[ii+zw]; // Pad right
    }
    memcpy(in,               in+zw+2,      (zw+2)*sizeof(float)); // Pad top
    memcpy(in+(zw+2)*(zh+1), in+(zw+2)*zh, (zw+2)*sizeof(float)); // Pad bottom
}

/// Replace each block of \a z x \a z pixels in \a in by its average.
/// Write the result in \a average.
/// Warning: image \a im is padded by 1 pixel each side. **TODO**
void projectionPU(const float*  in, int w, int h, int z,
                  float* temp, float* average) {
    float norm = 1.0f/float(z*z);
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int offset = w*j+i;
            float s=0.0;
            for(int y=j*z; y<(j+1)*z; y++)
                for(int x=i*z; x<(i+1)*z; x++){
                    int Zoffset = (z*w+2)*(y+1)+(x+1);
                    s += in[Zoffset];
                }
            s *= norm;
            temp[offset]=s;
        }
    zoom_duplication(temp,w,h, z, average);
}

void derivatives(const float a0[3], const float a1[3], const float a2[3],
                 float& grad, float& g, float& h) {
    const float sqrt12=M_SQRT1_2; // 1/sqrt(2)
    const float sqrt18=sqrt12/2;  // 1/sqrt(8)=1/(2 sqrt(2))
    const float c2=4.0*sqrt12+4.0;
    float c= a2[2] - a0[0];
    float d= a2[0] - a0[2];
    float ax= a1[2] - a1[0] + sqrt18*(c-d);
    float ay= a2[1] - a0[1] + sqrt18*(c+d);
    float az=ax*ay;
    ax *= ax;
    ay *= ay;
    grad = ax+ay;
    if(grad < ::minGrad) {
        h = g = ((a0[0]+a0[2]+a2[0]+a2[2])*sqrt12 +
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

/// Warning: image \a im is padded by 1 pixel each side.
void compute_derivatives(const float* im, int w, int h,
                         float* grad, float* evolxixi, float* evolnn) {
    // Inside
    const int w2 = w+2;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j=0; j<h; j++) {
        const float *a0=im+j*w2, *a1=a0+w2, *a2=a1+w2;
        float *gg=grad+j*w, *g=evolxixi+j*w, *h=evolnn+j*w;
        for(int i=0; i<w; i++) {
            derivatives(a0,a1,a2, *gg, *g, *h);
            a0++; a1++; a2++; gg++; g++; h++;
        }
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
    float* hrpad   = new float[iSizeImageZoom+2*(zw+zh+2)]; // Pad=1 each side

    const int n_iter = int(z*z/::dt+0.5f);
    std::cout << "Iterations: " << n_iter << "." << std::endl;

    // Generate dup, the 0-order spline interpolation of the input image
    zoom_duplication(lr,w,h, z, dup);
    // Copy it as initial zoomed image in padded image hrpad
    for(int i=0; i<zh; i++)
        memcpy(hrpad+(i+1)*(zw+2)+1, dup+i*zw, zw*sizeof(float));

    for(int k=0; k<n_iter; k++) {
        fill_padding(hrpad, zw, zh);
        compute_derivatives(hrpad,zw,zh, grad, uxixi, unn);
        projectionPU(hrpad,w,h, z, temp, average);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int y=0; y<zh; y++)
            for(int x=0; x<zw; x++) {
                int i=y*zw+x;
                float reac = average[i] - dup[i];
                float evolx = uxixi[i];
                float evoln = unn[i] / (1.0f + 0.01f * grad[i]);
                hrpad[(y+1)*(zw+2)+(x+1)] += ::dt*(evolx + evoln - reac);
                clamp(hrpad[(y+1)*(zw+2)+(x+1)]);
            }
    }

    // Remove padding
    for(int i=0; i<zh; i++)
        memcpy(hr+i*zw, hrpad+(i+1)*(zw+2)+1, zw*sizeof(float));

    delete [] temp;
    delete [] dup;
    delete [] uxixi;
    delete [] unn;
    delete [] grad;
    delete [] average;
    delete [] hrpad;
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
