#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j;
  int Sp; 
  double K_R, y_R; 

  Parameter_Table * Table = (Parameter_Table *)params;

  Sp  = Table->No_of_RESOURCES;
  y_R = Table->TOTAL_No_of_RESOURCES; /* Resource density (in number of resource units) */

  assert( Sp == 1 ); 

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  assert( Table->A == 0 );
  assert( Table->Lambda_C_0 ); 

  K_R = (double)Table->K_R; 

  for (j=0; j<Table->No_of_CELLS; j++) {

    A    = j*Table->LOCAL_STATE_VARIABLES + Table->A;
    RA   = j*Table->LOCAL_STATE_VARIABLES + Table->RA;

    dydt[A]  = -Table->Alpha_C_0 *y_R/K_R *y[A] + Table->Nu_C_0*y[RA];

    dydt[RA] = -Table->Nu_C_0*y[RA] + Table->Alpha_C_0 *y_R/K_R *y[A];
  }

  if(Table->No_of_CELLS > 1) { 
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