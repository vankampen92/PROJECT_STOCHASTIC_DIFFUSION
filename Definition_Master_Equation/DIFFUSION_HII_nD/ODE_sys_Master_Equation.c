#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function_ME (double t, const double y[], double dydt[], void *params)
{
  int i,j,k,s, No_of_CONFIGURATIONAL_STATES, D;
  int a_0, a_f, a_h;
  int a_fk, a_hk; 
  double AdUP, AdDW, OutUP, OutDW; 
  
  Parameter_Table * Table = (Parameter_Table *)params;
  
  Master_Equation * ME    = Table->MEq; 

  assert( Table->No_of_RESOURCES == ME->n_DIMENSION );
  assert( Table->No_of_CELLS == 1 ); 

  D     = Table->No_of_RESOURCES;
  a_0   = Table->TOTAL_No_of_CONSUMERS; 
  No_of_CONFIGURATIONAL_STATES = ME->No_of_CONFIGURATIONAL_STATES;

  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {

    AdUP = 0.0; 
    for(j=0; j < ME->Co[i]->anUp[D]; j++) {
      k = ME->Co[i]->anUp[j];
      
      for(s=0; s<Table->No_of_RESOURCES; s++) 
        if( (ME->Co[k]->n[s] - ME->Co[i]->n[s]) == 1 )
          AdUP += Table->Nu_Consumers[s] * (double)ME->Co[k]->n[s] * y[k];
    }

    AdDW = 0.0; 
    for(j=0; j < ME->Co[i]->anDw[D]; j++) {
      k = ME->Co[i]->anDw[j];

      a_hk = 0;
      for(s=0; s<Table->No_of_RESOURCES; s++) 
        a_hk += ME->Co[k]->n[s];
      a_fk = a_0 - a_hk;                       /* a_fk: Free consumers at the */
                                               /* k-th configuration state    */
      for(s=0; s<Table->No_of_RESOURCES; s++) 
        if( (ME->Co[k]->n[s] - ME->Co[i]->n[s]) == -1 )
          AdDW += Table->Theta_Consumers[s] * (double)a_fk * y[k];
    }

    a_h = 0;
    for(j=0; j<Table->No_of_RESOURCES; j++) 
      a_h += ME->Co[i]->n[j];
    a_f = a_0 - a_h;                           /* a_f: Free consumers at the */
                                               /* i-th configuration state   */
    OutUP = 0.0;
    for(j=0; j < ME->Co[i]->anUp[D]; j++) {
      k = ME->Co[i]->anUp[j];

      for(s=0; s<Table->No_of_RESOURCES; s++) 
        if( (ME->Co[k]->n[s] - ME->Co[i]->n[s]) == 1 )
          OutUP += Table->Theta_Consumers[s];
    }
    OutUP = OutUP * (double)a_f * y[i];

    OutDW = 0.0;
    for(j=0; j < ME->Co[i]->anDw[D]; j++) {
      k = ME->Co[i]->anDw[j];

      for(s=0; s<Table->No_of_RESOURCES; s++) 
        if( (ME->Co[k]->n[s] - ME->Co[i]->n[s]) == -1 )
            OutDW += Table->Nu_Consumers[s] * (double)ME->Co[i]->n[s];
    }
    OutDW = OutDW * y[i];   
    
    /* i-th Eq of the Master Equation for the i-th configuration state    */
	  dydt[i] = AdUP + AdDW - OutUP - OutDW;   /* dP(n_1, ..., n_S)/dt      */                                           
  }
  
  return GSL_SUCCESS;
}

int jacobian_ME (double t, const double y[], double *dfdy, double dfdt[],
		 void *params)
{
  int i;

  Parameter_Table * Table;
  gsl_matrix_view dfdy_mat;
  gsl_matrix * m;

  Table = (Parameter_Table *)params;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
  
  dfdy_mat = gsl_matrix_view_array(dfdy, K+1, K+1);
  
  m = &dfdy_mat.matrix;

  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[W], y[A]) */

  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[W]) */
  /* Optimally, this Function should be inline... */
  JACOBIAN_Matrix_ME(m, y, t, K, Table);
  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[W]) */

  return GSL_SUCCESS;
}

void JACOBIAN_Matrix_ME( gsl_matrix * m, const double *y, double t, int W_DUMMY,
			 Parameter_Table * Table)
{
  int i,k;
  double W_N, M,H, H_2, Q_Sigma_Deriv, Q_Recov_Deriv;
  double f; /* Infectious Humans */

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[W]) */
  /* First, setting entries to zero... */
  gsl_matrix_set_zero(m);
  
  // #include  <include.JAC_sys_ME.c>

  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[W]) */
}
