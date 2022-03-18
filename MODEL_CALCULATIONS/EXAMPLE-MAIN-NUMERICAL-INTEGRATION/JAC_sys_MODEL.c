/* 
	Jacobian of the System Defined 
*/
#include "MODEL.h"

int jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params)
{

  void JACOBIAN_Matrix(gsl_matrix *, const double *, double, int, Parameter_Model *);
  int i;

  Parameter_Model * Table;
  gsl_matrix_view dfdy_mat;
  gsl_matrix * m;

  Table = (Parameter_Model *)params;

  dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
  
  m = &dfdy_mat.matrix;

  JACOBIAN_Matrix(m, y, t, 3, Table);

  return GSL_SUCCESS;
}

void JACOBIAN_Matrix( gsl_matrix * m, const double *y, double t, int Dimension,
		      Parameter_Model * Table)
{
  double x; 
  
  /* First, setting entries to zero... */
  gsl_matrix_set_zero(m);

  double K_R;
  
  K_R = (double)Table->K_0;

  int R  = 0;
  int A  = 1;
  int RA = 2;
  
  /* F_R( y; parameters) */
  x= Table->Delta_R +Table->Beta_R -2.0*Table->Beta_R*y[R]/K_R -Table->Alpha_C *y[A]/K_R;
  gsl_matrix_set(m, R, 0, x);

  x= -Table->Alpha_C *y[R]/K_R;
  gsl_matrix_set(m, R, 1, x);

  x= 0.0;
  gsl_matrix_set(m, R, 2, x);

  /* F_A( y; parameters) */
  x= -Table->Alpha_C *y[A]/K_R;
  gsl_matrix_set(m, A, 0, x);

  x= -Table->Delta_C -Table->Alpha_C *y[R]/K_R ;
  gsl_matrix_set(m, A, 1, x);

  x= Table->Nu_C + Table->Beta_C;
  gsl_matrix_set(m, A, 2, x);

  /* F_RA( y; parameters) */
  x= Table->Alpha_C *y[A]/K_R;
  gsl_matrix_set(m, RA, 0, x);

  x= Table->Alpha_C *y[R]/K_R ;  
  gsl_matrix_set(m, RA, 1, x);

  x= -Table->Nu_C -Table->Delta_C;
  gsl_matrix_set(m, RA, 2, x);

  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[W]) */
}
