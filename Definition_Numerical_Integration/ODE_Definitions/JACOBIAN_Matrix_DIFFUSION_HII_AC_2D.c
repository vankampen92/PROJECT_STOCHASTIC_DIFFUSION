/*
   JACOBIAN_matrix_CASES_1.c is a companion file of JAC_sys_CASES_1.c
   and it is characterized by the following features:

   Jacobian corresponding to the ODE system implemented in ODE_sys_CASES_1.c:

   1. The recovery of CASES from clinical treatment. In ODE_sys_CASES_1.c, recovered individuals
   through treatment do no develop immunity at all, but reenter the susceptible class.

   2. Different possibilites to implement the dependence of the force of infection,
   $\beta$, on the mosquito infectious population, y[W], through a #define BETA_DERIVATIVE

   The Jacobian does not have any explicity dependence on time.

   Notice that JACOBIAN_matrix_CASES_1.c could be compiled
   independently from JAC_sys_CASES_1.c, but not the other
   way around.
*/
#include <MODEL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#define BETA_DERIVATIVE (Beta_Derivative_1(P,y[W]))

void JACOBIAN_Matrix( gsl_matrix * m, const double *y, double t, int W_DUMMY,
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
  
#include  <include.JAC_sys_DIFFUSION_HII_AC_2D.c>

  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[W]) */
}
