#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* This is the n Dim Holling Type II that corrsponds to a consumer feeding 
   on n different resources (No_of_RESOURCES)
*/

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j;
  int Sp; 
  double K_R, y_R, A_0; 

  Parameter_Table * Table = (Parameter_Table *)params;

  Sp  = Table->No_of_RESOURCES;

  assert( Table->Lambda_C_0 == 0.0 ); /* No addition of extra consumer through immmigration */
  assert( Sp == Table->LOCAL_STATE_VARIABLES );

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  assert( Table->A == 0 );

  A_0 = (double)Table->TOTAL_No_of_CONSUMERS; 
  K_R = (double)Table->K_R; 

  for (j=0; j<Table->No_of_CELLS; j++) {

    y_R = 0;
    for (i = 0; i < Table->LOCAL_STATE_VARIABLES; i++ ) {
      n    = j*Table->LOCAL_STATE_VARIABLES + i;
      y_R += y[n];   /*  Tota number of handling predators */
    }

    for (i = 0; i < Table->No_of_RESOURCES; i++ ) {
      n    = j*Table->LOCAL_STATE_VARIABLES + i;
      dydt[n]  = Table->Theta_Consumers[i] * (A_0 - y_R) - Table->Nu_Consumers[i]*y[n];
    }
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
