#include "MODEL.h"

int function (double t, const double y[], double dydt[], void * params)
{
  
  Parameter_Model * Table = (Parameter_Model *)params;

  int R    = 0; /* Resource */ 
  int A    = 1; /* Free Predators */ 
  int RA   = 2; /* Handling Predators */

  double K_R = Table->K_0; 

  dydt[R] = -Table->Delta_R *y[R] +Table->Beta_R *(K_R-y[R])/K_R *y[R] -Table->Alpha_C *y[R]/K_R *y[A];
    
  dydt[A] = -Table->Delta_C *y[A] +Table->Nu_C *y[RA] +Table->Beta_C*y[RA] -Table->Alpha_C *y[R]/K_R *y[A];

  dydt[RA] = -Table->Delta_C*y[RA] +Table->Alpha_C *y[R]/K_R *y[A] -Table->Nu_C*y[RA];

  return GSL_SUCCESS;
}


