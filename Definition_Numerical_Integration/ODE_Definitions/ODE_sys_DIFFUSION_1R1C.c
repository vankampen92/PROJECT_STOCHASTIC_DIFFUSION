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

  assert( Sp == 1 ); 

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; 

  for (j=0; j<Table->No_of_CELLS; j++) {

    R   = j*Table->LOCAL_STATE_VARIABLES + Table->R;
    A   = j*Table->LOCAL_STATE_VARIABLES + Table->A;
    RA  = j*Table->LOCAL_STATE_VARIABLES + Table->RA;
    ARA = j*Table->LOCAL_STATE_VARIABLES + Table->ARA;

    dydt[R] = Table->Lambda_R_0 *(K_R-y[R]) +Table->Beta_R *(K_R-y[R])/K_R *y[R] -Table->Delta_R_0 *y[R] -Table->Alpha_C_0 *y[R]/K_R *y[A];
    
    dydt[A] = Table->Lambda_C_0 -Table->Delta_C_0 *y[A] +2.0*Table->Nu_C_0 *y[RA] -Table->Alpha_C_0 *y[R]/K_R *y[A] -Table->Chi_C_0 *y[RA]/K_R *y[A] + Table->Eta_C_0 *y[ARA];

    dydt[RA] = Table->Alpha_C_0 *y[R]/K_R *y[A] -Table->Nu_C_0*y[RA] -Table->Chi_C_0 *y[RA]/K_R *y[A] +Table->Eta_C_0 *y[ARA];

    dydt[ARA] = Table->Chi_C_0 *y[RA]/K_R *y[A] -Table->Eta_C_0 *y[ARA]; 
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
