#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>


int function_ME (double t, const double y[], double dydt[], void *params)
{
  int i, No_of_CONFIGURATIONAL_STATES, n, m, a_0;
  int n_R; 
  int Sp; 
  double K_R, Theta, Nu, Chi, Eta; 
  int    i_1, i_2, i_3, i_4;
  double A_1, A_2, A_3, A_4, A_0; 
  
  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;
  assert( Sp == 1 );
  assert( Table->No_of_CELLS == 1 ); 

  K_R   = (double)Table->K_R;
  n_R   = Table->TOTAL_No_of_RESOURCES;   /* which determines resource density */
  a_0   = Table->TOTAL_No_of_CONSUMERS;   /* and the Theta parameter           */
  
  No_of_CONFIGURATIONAL_STATES = Table->MEq->No_of_CONFIGURATIONAL_STATES;

  Theta = Table->Alpha_C_0 * (double)n_R/K_R;  /* -H9  [Alpha_C_0] */
  Nu    = Table->Nu_C_0;                       /* -H10 [Nu_C_0]    */
  Chi   = Table->Chi_C_0;                      /* -H11             */                      
  Eta   = Table->Eta_C_0;                      /* -H12             */
  
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {

    i_to_nm_Map(Table, i, &n, &m);             /* n: A  (free predators)     */ 
                                               /* m: RA (handling predators) */

    nm_to_i_Map(Table, &i_1, n+1, m-1);
    nm_to_i_Map(Table, &i_2, n-1, m+1);
    nm_to_i_Map(Table, &i_3, n+1, m+1);
    nm_to_i_Map(Table, &i_4, n-1, m-1);
   
    A_1= Theta * (double)(n+1) * y[i_1];
    A_2= Nu *  (double)(m+1) * y[i_2];
    A_3= Chi * (double)(m+1)/K_R * (double)(n+1) * y[i_3];
    A_4= 0.5*Eta*(double)(a_0-n-m+2) * y[i_4];
    A_0= (Theta*(double)n+Nu*(double)m+Chi*(double)m/K_R*(double)n+0.5*Eta*(double)(a_0-n-m))*y[i];
    
    if( n == 0 ) { 
      if( m == 0 ) {
	dydt[i] = -A_0 + A_3;                    /* dP(n=0, m=0)/dt  */
      }
      else { 
	dydt[i] =  A_1 + A_3 - A_0;               /* dP(n=0, m)/dt    */
      }
    }
    else{
      if( m == 0 ) {
	dydt[i] = A_2 + A_3 - A_0;               /* dP(n, m=0)/dt    */
      }
      else {
	dydt[i] = A_1 + A_2 + A_3 + A_4 - A_0;   /* dP(n, m)/dt      */
	                                         /* n and m not zero */
      }
    }    
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
