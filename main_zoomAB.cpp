#include "cmdLine.h"
#include "io_png.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

const float sqrt2=M_SQRT1_2;
const float sqrt22=sqrt2/2;
const float seuil=0.001;

/// Replace each block of iSizeZoom x iSizeZoom pixels in pfZoom by its average.
/// Write the result in moyenne.
void ProjectionPU(const float *pfZoom, int dx, int dy, int iSizeZoom,
                  float *temp, float *moyenne) {
    for(int j=0; j<dy; j++)
        for(int i=0; i<dx; i++) {
            int offset = dx*j+i;
            float s=0.0;
            for(int x=i*iSizeZoom; x<(i+1)*iSizeZoom; x++)
                for(int y=j*iSizeZoom; y<(j+1)*iSizeZoom; y++)	{
                    int Zoffset = iSizeZoom*dx*y+x;
                    s += pfZoom[Zoffset];
                }
            s /= (float)(iSizeZoom*iSizeZoom);
            temp[offset]=s;
        }

    for(int j=0; j<dy; j++)
        for(int i=0; i<dx; i++) {
            int offset = dx*j+i;
            for(int x=i*iSizeZoom; x<(i+1)*iSizeZoom; x++)
                for(int y=j*iSizeZoom; y<(j+1)*iSizeZoom; y++) {
                    int Zoffset = (iSizeZoom*dx*y)+x;
                    moyenne[Zoffset]= temp[offset];
                }
        }
}

void Calcul(float *pict, int dx, int dy, float *evolxixi, float *evolnn, float *grad)
{
    int i,j,k;
    float c1,d1,l0,l1,l2,l3,l4,li,ax,ay,az,c2;
    int i1,i2,i3,i4,i5,i6;
    float *ai1,*ai2,*ai3,*ai4,*ai5,*ai6,*ai7,*ai8,*ai9,*gi5,*hi5,*ggi5;
    float *a,*g,*h,*gg,curv;
    int nx,ny,taille;
    /*-----------Input analysis --------------------*/
    nx=dx-1;
    ny=dy-1;
    taille = dx*dy;

    g=evolxixi;
    h=evolnn;
    gg=grad;
    a=pict;
    /*---------------------- PROCESS ---------------------------------------------*/
    /* ai1 .. ai9 are some pointers on our image pict */
    /*ai1 ai2 ai3 */
    /*ai4 ai5 ai6 */
    /*ai7 ai8 ai9 */

    c2=0.707106781186547129*4.0+4.0;
    ai1=a;
    ai2=a+1;
    ai3=a+2;
    ai4=a+(nx+1);
    ai5=ai4+1;
    gi5=g+(nx+2);
    hi5=h+(nx+2);
    ggi5=gg+(nx+2);
    ai6=ai5+1;
    ai7=ai4+(nx+1);
    ai8=ai7+1;
    ai9=ai8+1;
    curv = 0.0;
    /* Calculs principaux */   

    for(j=1;j<ny;j++){
        for(i=1;i<nx;i++){
        
            c1= *ai9- *ai1;
            d1= *ai7- *ai3; 
            ax= *ai6- *ai4 + sqrt22*(c1-d1);  
            ay= *ai8- *ai2 + sqrt22*(c1+d1);
            az=ax*ay;
            ax *= ax;
            ay *= ay;       	
            *ggi5= ax+ay;
            if (*ggi5<seuil) {
                *hi5=*gi5=(sqrt2*( *ai1+ *ai3+ *ai7+ *ai9)+( *ai2+ *ai4+ *ai6+ *ai8) - c2*( *ai5))/c2;
	   
            }
            else { 
                li=1.0/(*ggi5);
                az*=li;
                ax*=li;
                ay*=li;  
                li=az*az;
                l0=-2.0+4.0*li;
                l1=ay*(ay-ax);
                l2=ax*(ax-ay);
                l3=li-.5*az;
                l4=l3+az;
                *gi5=  l0*( *ai5)+l1*( *ai4+ *ai6)+l2*( *ai2+ *ai8)+l3*( *ai9+ *ai1)+l4*( *ai3+ *ai7);
	   
 
                *hi5=l0*( *ai5)+l2*( *ai4+ *ai6)+l1*( *ai2+ *ai8)+l4*( *ai9+ *ai1)+l3*( *ai3+ *ai7); 
           
            }
            ai1=ai2; ai2=ai3; ai3++; ai4=ai5; ai5=ai6; ai6++; ai7=ai8; ai8=ai9; ai9++; gi5++; hi5++; ggi5++;
        }
        ai1=ai2; ai2=ai3; ai3++; ai4=ai5; ai5=ai6; ai6++; ai7=ai8; ai8=ai9; ai9++; gi5++; hi5++; ggi5++;
        ai1=ai2; ai2=ai3; ai3++; ai4=ai5; ai5=ai6; ai6++; ai7=ai8; ai8=ai9; ai9++; gi5++; hi5++; ggi5++;	
    }
    /* Les bords ... */
    i1=1; i2=nx+2; i3=ny*(nx+1)+1; i4=i3-nx-1; 
    for(k=1;k<nx;k++){  
        c1=a[i2+1]-a[i1];
        d1=a[i2-1]-a[i1];
        ax=((a[i1+1]-a[i1-1])+sqrt22*(c1-d1));
        ay=(2.0*(a[i2]-a[i1])+sqrt22*(c1+d1));
        az=ax*ay;
        ax *= ax;
        ay *= ay;
        gg[i1]=ax+ay;
        if (gg[i1]<seuil) 
            { 
                h[i1]= g[i1]=((a[i1-1]+a[i1+1]+a[i2])*2.0+a[i2-1]+a[i2+1]-8.0*a[i1])/8.0; 
            }
        else {
            li=1.0/gg[i1];
            az*=li;
            ax*=li;
            ay*=li; 
            li=az*az;
            l0=-2.0+4.0*li;
            l2=ax*(ax-ay);
            l1=ay*(ay-ax);
            l3=li;	
            g[i1]=l0*a[i1]+l1*(a[i1-1]+a[i1+1])+l2*(a[i2]+a[i2])+(l3+l3)*(a[i2-1]+a[i2+1]);
            h[i1]=l0*a[i1]+l2*(a[i1-1]+a[i1+1])+l1*(a[i2]+a[i2])+(l3+l3)*(a[i2-1]+a[i2+1]);
        }

        c1=a[i3]-a[i4-1];
        d1=a[i3]-a[i4+1];
        ax=((a[i3+1]-a[i3-1])+sqrt22*(c1-d1));
        ay=(2.0*(a[i3]-a[i4])+sqrt22*(c1+d1));
        az=ax*ay;
        ax *= ax;
        ay *= ay;
        gg[i3]=ax+ay;
        if (gg[i3]<seuil) {
            h[i3]=g[i3]=((a[i3-1]+a[i3+1]+a[i4])*2.0+a[i4-1]+a[i4+1]-8.0*a[i3])/8.0;	
        }
        else {
            li=1.0/gg[i3];
            az*=li;
            ax*=li;
            ay*=li; 
            li=az*az;
            l0=-2.0+4.0*li;
            l2=ax*(ax-ay);
            l1=ay*(ay-ax);
            l3=li;
            g[i3]=l0*a[i3]+l1*(a[i3-1]+a[i3+1])+l2*(a[i4]+a[i4])+(l3+l3)*(a[i4-1]+a[i4+1]);
            h[i3]=l0*a[i3]+l2*(a[i3-1]+a[i3+1])+l1*(a[i4]+a[i4])+(l3+l3)*(a[i4-1]+a[i4+1]);
	
        }	
        i1++; i2++; i3++; i4++;
    }
    i1=0; i2=nx+1; i3=i2+i2; i4=nx; i5=nx+nx+1; i6=i5+nx+1;
    for(k=1;k<ny;k++){ 
        c1=a[i3+1]-a[i2];
        d1=a[i2]-a[i1+1];
        ax=(2.0*(a[i2+1]-a[i2])+sqrt22*(c1-d1));
        ay=((a[i3]-a[i1])+sqrt22*(c1+d1));
        az=ax*ay;
        ax *= ax;
        ay *= ay;
        gg[i2]=ax+ay;   
        if (gg[i2]<seuil) { 
            h[i2]= g[i2]=((a[i1]+a[i3]+a[i2+1])*2.0+a[i3+1]+a[i1+1]-8.0*a[i2])/8.0;	
        }
        else {
            li=1.0/gg[i2];
            az*=li;
            ax*=li;
            ay*=li;
            li=az*az;
            l0=-2.0+4.0*li;
            l2=ax*(ax-ay);
            l1=ay*(ay-ax);
            l3=li;	
            g[i2]=l0*a[i2]+l2*(a[i3]+a[i1])+l1*(a[i2+1]+a[i2+1])+(l3+l3)*(a[i1+1]+a[i3+1]);
            h[i2]=l0*a[i2]+l1*(a[i3]+a[i1])+l2*(a[i2+1]+a[i2+1])+(l3+l3)*(a[i1+1]+a[i3+1]);
        }	
        c1=a[i5]-a[i4-1];
        d1=a[i6-1]-a[i5];
        ax=(2.0*(a[i5]-a[i5-1])+sqrt22*(c1-d1));
        ay=((a[i6]-a[i4])+sqrt22*(c1+d1));
        az=ax*ay;
        ax *= ax;
        ay *= ay;
        gg[i5]=ax+ay; 
        if (gg[i5]<seuil) { 
            h[i5]= g[i5]=((a[i4]+a[i6]+a[i5-1])*2.0+a[i4-1]+a[i6-1]-8.0*a[i5])/8.0;
	
        }
        else {
            li=1.0/gg[i5];
            az*=li;
            ax*=li;
            ay*=li; 
            li=az*az;
            l0=-2.0+4.0*li;
            l2=ax*(ax-ay);
            l1=ay*(ay-ax);
            l3=li;	
            g[i5]=l0*a[i5]+l2*(a[i6]+a[i4])+l1*(a[i5-1]+a[i5-1])+(l3+l3)*(a[i4-1]+a[i6-1]);
            h[i5]=l0*a[i5]+l1*(a[i6]+a[i4])+l2*(a[i5-1]+a[i5-1])+(l3+l3)*(a[i4-1]+a[i6-1]);
        }   
        i1=i2; i2=i3; i3+=nx+1; i4=i5; i5=i6; i6+=nx+1;
    }
   	h[0]=g[0]=(a[1]+a[nx+2]+a[nx+1]-3.0*a[0])/6.0;        
   	h[nx]=g[nx]=(a[nx-1]+a[nx+nx]+a[nx+nx+1]-3.0*a[nx])/6.0;	
    i1=ny*(nx+1);
   	h[i1]=g[i1]=(a[i1-nx-1]+a[i1-nx]+a[i1+1]-3.0*a[i1])/6.0;      
    i1=i1+nx;
   	h[i1]=g[i1]=(a[i1-1]+a[i1-nx-1]+a[i1-nx-2]-3.0*a[i1])/6.0;
}

void zoomAB(const float* pfImage, int dx, int dy, float* pfZoom, int iSizeZoom) {
    const int zdx=iSizeZoom*dx, zdy=iSizeZoom*dy;
    const int iSizeImageZoom = zdx*zdy;
    float *pfImagetemp = new float [dx*dy];
    float *pfImageDup = new float [iSizeImageZoom];
    float *pfuxixi = new float [iSizeImageZoom];
    float *pfunn   = new float [iSizeImageZoom];
    float *pfgrad  = new float [iSizeImageZoom];
    float *pfmoyenne  = new float [iSizeImageZoom];

    float dt = (float)0.1;
	int n_iter= (int) ((double(iSizeZoom*iSizeZoom))/dt+0.5);

    std::cout << "Iterations: " << n_iter << std::endl;
	// Ici Dupliquer pfImage vers pfImageDup
	for (int i=0;i<dx;i++) 
		for (int j=0;j<dy;j++) {
            int offset = dx*j+i;
            for (int x=i*iSizeZoom; x<(i+1)*iSizeZoom; x++)
                for (int y=j*iSizeZoom; y<(j+1)*iSizeZoom; y++) {
                    int Zoffset = (zdx*y)+x;
                    pfImageDup[Zoffset]= pfImage[offset];
                }
        }
	// copier pfImageDup vers pfZoom
    memcpy(pfZoom, pfImageDup, iSizeImageZoom * sizeof(float));

    for(int k=0; k<n_iter; k++) {
        std::cout << k << ' ' << std::flush;
        Calcul(pfZoom,zdx,zdy,pfuxixi,pfunn,pfgrad);
        ProjectionPU(pfZoom,dx,dy,iSizeZoom,pfImagetemp,pfmoyenne);

		const float* uxi=pfuxixi;
        const float* un=pfunn;
        const float* grd=pfgrad;
        const float* moy= pfmoyenne;
        const float* uo = pfImageDup; 
		float* imz=pfZoom;
        std::cout << "(" << std::flush;
        for(int i=0; i<iSizeImageZoom; i++, imz++) {
            std::cout << i << ' ' << std::flush;
            float reac = *moy++ - *uo++;
            std::cout << reac << ' ' << std::flush; // DEBUG
            float evolx = *uxi++;
            std::cout << evolx << ' ' << std::flush; // DEBUG
            std::cout << *un << ' ' << std::flush; // DEBUG
            std::cout << *grd << ' ' << std::flush; // DEBUG
            float evoln = *un++ / (1.0 + 0.01 * *grd++);
            std::cout << evoln << ' ' << std::flush; // DEBUG
            *imz += (evolx + evoln - reac)*(dt);
            if(*imz > 255.0)
				*imz = 255;
            else if(*imz <0.0)
                *imz = 0;
            break; // DEBUG
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
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
