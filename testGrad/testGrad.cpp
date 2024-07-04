#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>

const int RADIUS=10; ///< Radius of disk
const int SIZE=2*RADIUS+10;   ///< Dimension of square image

#define straux(x) #x
#define str(x) straux(x)

// Choices of k00 for curvature computation. Select one and comment the others.
//#define K00 k00 = 0.5*d2 // Sapiro-Tannenbaum
//#define K00 k00 = 0.25*d2 // Alvarez
#define K00 k00 = 0.5*d2-gx*gx*gy*gy/d2 // Alvarez-Morel
//#define K00 k00 = (d2-fabs(gx*gy)+fmax(gx*gx,gy*gy))/4 // Cohignac et al.
//#define K00 k00 = (d2-fabs(gx*gy)+2*fmax(gx*gx,gy*gy))/6 // Monasse

/// Return Gaussian kernel of size 2*radius+1. \a radius is set to ceil(3*sigma)
/// with sigma a constant.
float* gaussKernel(int& radius) {
    const float SIGMA=1.6f;
    radius = (int)ceil(3*SIGMA);
    int sz = 2*radius+1;
    float sum=0, *ker = new float[sz];
    for(int i=0; i<=radius; i++) {
        ker[sz-i-1] = ker[i] = exp(-(i-radius)*(i-radius)/(2*SIGMA*SIGMA));
    }
    for(int i=0; i<sz; i++)
        sum += ker[i];
    for(int i=0; i<sz; i++)
        ker[i] /= sum;
    return ker;
}

/// x-convolution and 1/2 x-subsampling. The result is transposed.
void convx_sub2(const float* im, int w, int h, const float* ker, int radius,
                float*& out) {
    assert(radius<=2*w);
    out = new float[h*(w/2)];
    for(int y=0; y<h; y++) {
        const float* in=im+y*w;
        for(int x=0; x<w; x+=2) {
            float c=0;
            const float* k = ker;
            for(int i=-radius; i<=+radius; i++) {
                int j = x+i;
                if(j<0) j = -j-1;
                if(j>=w) j = 2*w-j-1;
                c += in[j] * *k++;
            }
            out[(x/2)*h+y] = c;
        }
    }
}

/// Compute gradient of image at (x,y) with stencils
/// (-C 0 C)               (-C -1 -C)
/// (-1 0 1) / 2(2C+1) and ( 0  0  0) / 2(2C+1)
/// (-C 0 C)               ( C  1  C)
void grad(const float* im, int x, int y, float& gx, float& gy, float C) {
    im += x+y*SIZE;
    gx = im[-1-SIZE]*(-C) + im[+1-SIZE]*C +
         im[-1]     *(-1) + im[1]      *1 +
         im[-1+SIZE]*(-C) + im[+1+SIZE]*C;
    gy = im[-1-SIZE]*(-C) + im[0-SIZE]*(-1) + im[+1-SIZE]*(-C) +
         im[-1+SIZE]*(+C) + im[0+SIZE]*(+1) + im[+1+SIZE]*(+C);
    float norm=2*(2*C+1);
    gx /= norm;
    gy /= norm;
}

float curv(const float* im, int x, int y, float gx, float gy) {
    double d2 = gx*gx+gy*gy;
    float K00; // See macro at top of file for choice of k00
    float k10 = 2*k00-gy*gy;
    float k01 = 2*k00-gx*gx;
    float k1_1 = -k00+(d2+gx*gy)/2;
    float k11  = -k00+(d2-gx*gy)/2;
    im += x+y*SIZE;
    return (-4*k00*im[0] + k10*(im[-1]+im[1]) + k01*(im[-SIZE]+im[SIZE])
            + k1_1*(im[-SIZE+1]+im[SIZE-1]) + k11*(im[-SIZE-1]+im[SIZE+1])) /d2;
}

/// Return average and standard deviation of signal.
void stats(const std::vector<float>& v, float& av, float& std) {
    std::vector<float>::const_iterator it, end=v.end();
    av=0;
    for(it=v.begin(); it!=end; ++it)
        av += *it;
    av /= v.size();
    std=0;
    for(it=v.begin(); it!=end; ++it)
        std += (*it-av)*(*it-av);
    std /= v.size();
    std = sqrt(std);
}

int main(int argc, char* argv[]) {
    float C=1/sqrt(2);
    if(argc>1)
        C = std::stof(argv[1]);

    // Generate image
    float* im = new float[2*SIZE * 2*SIZE];
    float* p=im;
    for(int i=-SIZE; i<SIZE; i++)
        for(int j=-SIZE; j<SIZE; j++)
            *p++ = hypot(i,j)<2*RADIUS? 0: 1;

    // Reduce /2 image
    int radius;
    float* ker = gaussKernel(radius);
    float* xconv;
    convx_sub2(im, 2*SIZE, 2*SIZE, ker, radius, xconv);
    delete [] im;
    convx_sub2(xconv, 2*SIZE, SIZE, ker, radius, im);
    delete [] xconv;
    delete [] ker;

    // Measure
    const int length=(int)ceil(2*M_PI*RADIUS);
    std::vector<float> modulus, angle, curvature;
    for(int i=0; i<length; i++) {
        float theta = i*2*M_PI/length;
        float dx=cos(theta), dy=sin(theta);
        float m = fmaxf(fabs(dx),fabs(dy));
        dx /= m;
        dy /= m;
        float g=0, alpha=0, c=0;
        int jopt = -1;
        for(int j=RADIUS-0; j<=RADIUS+0; j++) { // Test close to discontinuity
            int x = SIZE/2+(int)round(j*dx), y = SIZE/2+(int)round(j*dy);
            if(! (0<=x-1 && x+1<SIZE && 0<=y-1 && y+1<SIZE))
                break;
            float gx, gy;
            grad(im, x, y, gx, gy, C);
            float d = hypot(gx,gy);
            if(g<d) {
                jopt = j;
                g=d;
                alpha = m*dx * gx/d + m*dy * gy/d;
                c = curv(im, x, y, gx, gy);
            }
        }
        assert(g>0);
        std::cout << jopt << ' ';
        modulus.push_back(g);
        angle.push_back(acos(alpha));
        curvature.push_back(c);
    }
    std::cout << std::endl;
    std::cout << "C=" << C << ", " << str(K00) << std::endl;
    float m, std;
    stats(modulus, m, std);
    std::cout << "Modulus: " << m << " +/- " << std << std::endl;
    stats(angle, m, std);
    std::cout << "Angle:   " << m << " +/- " << std << std::endl;
    stats(curvature, m, std);
    std::cout << "Curv:    " << m << " +/- " << std << std::endl;
    delete [] im;
}
