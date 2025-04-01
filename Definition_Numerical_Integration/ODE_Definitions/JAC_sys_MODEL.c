/* 
	General File (valid for all models)
*/

#include <MODEL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#define BETA_DERIVATIVE (Beta_Derivative_1(P,y[W]))

int 
jacobian (double t, const double y[], double * dfdy, double dfdt[],
	        void *params)
{

  void JACOBIAN_Matrix(gsl_matrix *, const double *, double, int, Parameter_Table *);
  int i;

  Parameter_Table * Table;
  gsl_matrix_view dfdy_mat;
  gsl_matrix * m;

  Table = (Parameter_Table *)params;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
  
  dfdy_mat = gsl_matrix_view_array(dfdy, K+1, K+1);
  
  m = &dfdy_mat.matrix;

  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  /* Optimally, this Function should be inline... */
  JACOBIAN_Matrix(m, y, t, K, Table);
  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */

  return GSL_SUCCESS;
}

