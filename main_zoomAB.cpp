#include "cmdLine.h"
#include "io_png.h"
#include <iostream>
#include <cmath>
#include <cstring>

const float sqrt2=M_SQRT1_2;
const float sqrt22=sqrt2/2;
const float minGrad=0.001f; // Minimum gradient norm for applying anisotropy
const float dt = 0.1f; // Time step for diffusion

/// Replace each block of iSizeZoom x iSizeZoom pixels in pfZoom by its average.
/// Write the result in moyenne.
void projectionPU(const float *pfZoom, int w, int h, int iSizeZoom,
                  float *temp, float *moyenne) {
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int offset = w*j+i;
            float s=0.0;
            for(int x=i*iSizeZoom; x<(i+1)*iSizeZoom; x++)
                for(int y=j*iSizeZoom; y<(j+1)*iSizeZoom; y++) {
                    int Zoffset = iSizeZoom*w*y+x;
                    s += pfZoom[Zoffset];
                }
            s /= (float)(iSizeZoom*iSizeZoom);
            temp[offset]=s;
        }

    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int offset = w*j+i;
            for(int x=i*iSizeZoom; x<(i+1)*iSizeZoom; x++)
                for(int y=j*iSizeZoom; y<(j+1)*iSizeZoom; y++) {
                    int Zoffset = (iSizeZoom*w*y)+x;
                    moyenne[Zoffset]= temp[offset];
                }
        }
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
    const float *a0=im, *a1=a0+w, *a2=a1+w;     
    float *gg=grad+w+1, *g=evolxixi+w+1, *h=evolnn+w+1;
    for(int j=1; j+1<dy; j++){
        for(int i=1; i+1<w; i++){
            derivatives(a0,a1,a2, *gg, *g, *h);
            a0++; a1++; a2++; gg++; g++; h++;
        }
        a0+=2; a1+=2; a2+=2; gg+=2; g+=2; h+=2;
    }
    // Top boundary
    a0=im; a1=im+w;
    gg=grad; g=evolxixi; h=evolnn;
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

void zoomAB(const float* pfImage, int w, int h, float* pfZoom, int iSizeZoom) {
    const int zw=iSizeZoom*w, zh=iSizeZoom*h;
    const int iSizeImageZoom = zw*zh;
    float *pfImagetemp = new float[w*h];
    float *pfImageDup = new float[iSizeImageZoom];
    float *pfuxixi =    new float[iSizeImageZoom];
    float *pfunn   =    new float[iSizeImageZoom];
    float *pfgrad  =    new float[iSizeImageZoom];
    float *pfmoyenne  = new float[iSizeImageZoom];

	const int n_iter = int(iSizeZoom*iSizeZoom/::dt+0.5f);

    std::cout << "Iterations: " << n_iter << std::endl;
	// Ici Dupliquer pfImage vers pfImageDup
	for (int i=0;i<w;i++)
		for (int j=0;j<h;j++) {
            int offset = w*j+i;
            for (int x=i*iSizeZoom; x<(i+1)*iSizeZoom; x++)
                for (int y=j*iSizeZoom; y<(j+1)*iSizeZoom; y++) {
                    int Zoffset = (zw*y)+x;
                    pfImageDup[Zoffset]= pfImage[offset];
                }
        }
	// copier pfImageDup vers pfZoom
    memcpy(pfZoom, pfImageDup, iSizeImageZoom * sizeof(float));

    for(int k=0; k<n_iter; k++) {
        std::cout << k << ' ' << std::flush;
        compute_derivatives(pfZoom,zw,zh, pfgrad, pfuxixi, pfunn);
        projectionPU(pfZoom,w,h, iSizeZoom, pfImagetemp, pfmoyenne);

		const float* uxi=pfuxixi;
        const float* un=pfunn;
        const float* grd=pfgrad;
        const float* moy= pfmoyenne;
        const float* uo = pfImageDup; 
		float* imz=pfZoom;
        for(int i=0; i<iSizeImageZoom; i++, imz++) {
            float reac = *moy++ - *uo++;
            float evolx = *uxi++;
            float evoln = *un++ / (1.0f + 0.01f * *grd++);
            *imz += ::dt*(evolx + evoln - reac);
            if(*imz > 255.0)
				*imz = 255;
            else if(*imz <0.0)
                *imz = 0;
        }
    }
    std::cout << std::endl;
    delete [] pfImagetemp;
    delete [] pfImageDup;
    delete [] pfuxixi;
    delete [] pfunn;
    delete [] pfgrad;
    delete [] pfmoyenne;
}

int main(int argc, char* argv[]) {
    CmdLine cmd;
    unsigned int z=2; // Max number of pixels to discard
    cmd.add( make_option('z',z).doc("Integer zoom factor") );

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

    size_t w, h;
    float* data = io_png_read_f32_gray(argv[1], &w, &h);
    if(! data) {
        std::cerr << "Error loading image " << argv[1] << std::endl;
        return 1;
    }

    float* out = new float[1*z*z*w*h];
    zoomAB(data, (int)w, (int)h, out, (int)z);
    if(io_png_write_f32(argv[2], out, z*w, z*h, 1) != 0) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }

    free(data);
    delete [] out;
    return 0;
}
