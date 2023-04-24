/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"
#include <stdlib.h>
#include <math.h>

void tridag_Function( double (*a)(int, int), 
		      double (*b)(int, int),
		      double (*c)(int, int), 
		      double (*r)(int, int), 
		      double * u,
		      unsigned long n )
{
  unsigned long j;
  double bet;
  double * gam = (double *)malloc( n * sizeof(double));
	
  if ((*b)(0,n) == 0.0) nrerror("Error 1 in tridag");

  bet=(*b)(0,n);
  u[0]= (*r)(0,n) / bet;

	for (j=1;j < n;j++) {
	  gam[j] = (*c)(j-1,n) / bet;
	  bet= (*b)(j,n) - (*a)(j,n)*gam[j];
	  if (bet == 0.0)	nrerror("Error 2 in tridag");
	  u[j] = ((*r)(j,n)- (*a)(j,n) * u[j-1])/bet;
	}
	
	for (j=(n-2); j>=0 ;j--)
	  u[j] -= gam[j+1]*u[j+1];

  free(gam);
}
#undef NRANSI
