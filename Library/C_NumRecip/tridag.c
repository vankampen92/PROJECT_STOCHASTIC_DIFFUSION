/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"
#include <stdlib.h>
#include <math.h>

void tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n)
{
  /* Single precision function.
     Fortran inspired treatment for indexing vectors:
     a[1],...,a[N]
  */
	unsigned long j;
	float bet,*gam;

	gam=vector(1,n);
	if (b[1] == 0.0) nrerror("Error 1 in tridag");
	u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0)	nrerror("Error 2 in tridag");
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
	free_vector(gam,1,n);
}

void tridag_double(double a[], double b[], double c[], double r[], double u[],
		   unsigned long n)
{
  /* Double precision function.
     C inspired treatment for indexing vectors:
     a[0],...,a[N-1]
  */
  /*  Solution of a tridag linear system: 

    |   b[0]   c[0]   0    0                 0        | | u[0]   |     | r[0]   | 
    |   a[1]   b[1]  c[1]  0                 0        | |        |     |        |
    |    0     a[2]  b[2] c[2]               0        | |        |     |        |
    |                                        0        | |        |  =  |        |
    |                                        0        | |        |     |        |
    |                                      c[n-2]     | |        |     |        |
    |    0      0     0    0       a[n-1]  b[n-1]     | | u[n-1] |     | r[n-1] |

    where imputs are a[], b[], c[], defining the matrix,
    and r[] defining the indepenent vector. The output is given in u[]:
  */
	int j;
	double bet;
	double * gam = (double *)malloc( n * sizeof(double));

	if (b[0] == 0.0) nrerror("Error 1 in tridag");
	u[0]=r[0]/(bet=b[0]);
	for (j=1;j<n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0)	nrerror("Error 2 in tridag");
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}

	for (j=(n-2);j>=0;j--) {
	  u[j] -= gam[j+1]*u[j+1];
	}

	free(gam);
}

#undef NRANSI
