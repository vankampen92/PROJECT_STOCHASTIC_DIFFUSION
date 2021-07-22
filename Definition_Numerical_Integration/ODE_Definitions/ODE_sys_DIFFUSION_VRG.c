#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j;
  int Sp; 
  double K_R; 

  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;

  assert( Sp == 1 );         // Only one vegetation type

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; 

  for (j=0; j<Table->No_of_CELLS; j++) {

    V   = j*Table->LOCAL_STATE_VARIABLES + Table->V;
    R   = j*Table->LOCAL_STATE_VARIABLES + Table->R;
    G   = j*Table->LOCAL_STATE_VARIABLES + Table->G;

    dydt[V] = Table->Lambda_R_0 *(K_R-y[V]-y[G]) + (Table->Beta_R + Table->p_1*y[G]) *(K_R-y[V]-y[G])/K_R *y[V] - Table->Delta_R_0 *y[V] - Table->Alpha_C_0 *y[V]/K_R *y[R];
    
    dydt[R] = Table->Lambda_C_0*K_R         + Table->Theta_C*Table->Alpha_C_0 *y[V]/K_R *y[R]         - Table->Delta_C_0 *y[R];

    dydt[G] = Table->Lambda_C_1 *(K_R-y[V]-y[G]) + Table->Beta_C *(K_R-y[V]-y[G])/K_R* y[G]           - Table->Delta_C_1 *y[G];

  }

  if( Table->No_of_CELLS > 1) {   
    n= 0; 
    for (j=0; j<Table->No_of_CELLS; j++) { 
   
      for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++) { 
      dydt[n] += In_Mu(Table, n, i, j, y) - Out_Mu_Per_Capita(Table, i, j) * y[n];;
      n++;
      }
    }
  }
  
  return GSL_SUCCESS;
}
