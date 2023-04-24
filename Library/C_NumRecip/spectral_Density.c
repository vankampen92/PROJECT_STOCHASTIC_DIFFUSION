#include <math.h>
#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))       /* Bartlett */

void spectral_Density(float data[], float p[], int m, int k, int ovrlap)
{
  /* 
     Returns data in p[] as the data's power (mean square amplitude) at 
     frequency $(j-1)/(2*m)$ cycles per grid point, for j = 1 to m+1.
     (f_c = 1/2). The number of data in w[] must be greater than 
     (2*k+1)*m if ovrlap = 1 and (4*k*m) if ovrlap = 0. The number 
     of segments of the data is 2*m in both cases. The routine calls
     four1 k times, each call with 2 partitions of 2*m points 
  */
  void four1(float data[], unsigned long nn, int isign);
  int no, mm,m44,m43,m4,kk,joffn,joff,j2,j;
  float w,facp,facm,*w1,*w2,sumw=0.0,den=0.0;

	mm=m+m;
	m43=(m4=mm+mm)+3;
	m44=m43+1;
	w1=vector(1,m4);
	w2=vector(1,m);
	facm=m;
	facp=1.0/m;
	for (j=1;j<=mm;j++) sumw += SQR(WINDOW(j,facm,facp));
	for (j=1;j<=m+1;j++) p[j]=0.0;
	no = 0;
	if (ovrlap)
	  for (j=1;j<=m;j++){
	    no++;
	    w2[j] = data[j];
	  }  
	/* Loop over data segments in groups of two intertwined segments */
	for (kk=1;kk<=k;kk++) {
		for (joff = -1;joff<=0;joff++) {
		        if (ovrlap) {
				for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];
				for (j=1;j<=m;j++) w2[j] = data[no+j];
				no+=m;
				joffn=joff+mm;
				for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];
			} else {
			        for (j=joff+2;j<=m4;j+=2) w1[j] = data[++no];
			}
		}
		/* Apply the window to the data */
		for (j=1;j<=mm;j++) {
			j2=j+j;
			w=WINDOW(j,facm,facp);
			w1[j2] *= w;
			w1[j2-1] *= w;
		}
		four1(w1,mm,1);
		p[1] += (SQR(w1[1])+SQR(w1[2]));
		for (j=2;j<=m;j++) {
			j2=j+j;
			p[j] += (SQR(w1[j2])+SQR(w1[j2-1])
				+SQR(w1[m44-j2])+SQR(w1[m43-j2]));
		}
		p[m+1] += (SQR(w1[m+1])+SQR(w1[m+2]));
		den += sumw;
	}
	den *= m4;           /* Correct normalization factor!!! */ 
	for (j=1;j<=m+1;j++) p[j] /= den;  
	free_vector(w2,1,m);
	free_vector(w1,1,m4);
}
#undef WINDOW
#undef NRANSI
