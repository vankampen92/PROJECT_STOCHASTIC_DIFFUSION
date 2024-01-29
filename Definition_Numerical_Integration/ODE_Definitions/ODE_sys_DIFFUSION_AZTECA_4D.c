#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j;
  int Sp; 

  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;

  assert( Sp == 1 ); 

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  /* K_W: Total Carrying Capacity (Workers) per local cell, 
     where K_R is the max No of Workers per Nest */
  double K_W = (double)Table->K_R * (double)Table->Lambda_C_1;
  /* K_Q: Total Max No of Nests (per local patch) */ 
  double K_Q = (double)Table->Lambda_C_1; 

  for (j=0; j<Table->No_of_CELLS; j++) {

    W  = j*Table->LOCAL_STATE_VARIABLES + Table->W;
    Q   = j*Table->LOCAL_STATE_VARIABLES + Table->Q;
    F   = j*Table->LOCAL_STATE_VARIABLES + Table->F;
    WF  = j*Table->LOCAL_STATE_VARIABLES + Table->WF;

    /* Lambda_R_0 represents an exteral passive
       source of workder arriving in the local
       patches.
       Lambda_C_0 represents an exteral passive
       source of flies arriving in the local
       patches.
    */
    dydt[W]  = Table->Lambda_R_0 -Table->Delta_R_0 *y[W] +Table->Beta_R *(K_W - y[W])/K_W * y[Q] -Table->Alpha_C_0*y[W]/K_W *y[F];             
    
    dydt[Q]  = Table->Eta_R *(K_Q-y[Q])/K_Q *y[W] -Table->Delta_R_1*y[Q]; 
    
    dydt[F]  = Table->Lambda_C_0 + Table->Nu_C_0 *y[WF] -Table->Delta_C_0 *y[F] ;

    dydt[WF] = Table->Alpha_C_0*y[W]/K_W *y[F] -Table->Delta_C_1*y[WF] -Table->Nu_C_0*y[WF] ;

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
