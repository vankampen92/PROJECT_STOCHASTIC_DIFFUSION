#include "../Include/MODEL.h"

#define BIOLOGICAL_ZERO_POPULATION (-5.0e-1)
#define GSL_STEP_RK4 (0.01)        /* Integration step: 1/100 fraction of a day */
                                   /* for non-adaptive step integration     */
#if defined GSL_RK4
#define NON_ADAPTIVE_STEP_INTEGRATION
#endif 

#define GSL_STEP_MINIM 1.0e-10 /* 1/1.e+20 fraction of a day */


#define GSL_ODE_ERR_0 (1.0e-5)  // Original value from GSL : 1.0e-5
#define GSL_ODE_ERR_1 (1.0e-6)  // Original value form GSL : 1.0e-6

int master_equation_driver( Parameter_Table * Table, 
			    int j, double * Time_Current )
{
  int status;
#if defined GSL_BSIMP
  const gsl_odeiv_step_type * T = gsl_odeiv_step_bsimp; /* This 'bsimp' uses the Jacobian */
#else
#if defined GSL_RK4
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
#else
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
#endif
#endif

  int MODEL_STATE_VARIABLES = Table->MEq->No_of_CONFIGURATIONAL_STATES;
  /* Total No of Discrete Values in Probability Distribution, i. e., 
     the Total No of Equations in the Master Equaiton System. 
     For systems that have naturally an infinit number of configurational states, 
     the Total No of Equations in the Master Equation, this is, the total number 
     of configurational states should be caped at some large enough value. In 
     this case, the ME integration is only an approximation of the true 
     infinit system. 
  */
  
  gsl_odeiv_step * s 
         = gsl_odeiv_step_alloc (T, MODEL_STATE_VARIABLES);
  gsl_odeiv_control * c 
         = gsl_odeiv_control_y_new (GSL_ODE_ERR_0, GSL_ODE_ERR_1);
  gsl_odeiv_evolve * e 
         = gsl_odeiv_evolve_alloc (MODEL_STATE_VARIABLES);
     
  gsl_odeiv_system sys = {function_ME, jacobian_ME, MODEL_STATE_VARIABLES, Table};
     
  double * t = Time_Current;
  double t1  = Table->T->Time_Vector[j];

  double h = 1e-6;
  
  double * y  = Table->MEq->Probability_Distribution;
    
  while (   (*t) < t1   )
    {
      status = gsl_odeiv_evolve_apply (e, c, s,
				       &sys, 
				       t, t1,
				       &h, y);
     
      if (status != GSL_SUCCESS) break;
    }
     
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  
  
  return status;
}

void Normalization_Master_Equation(Parameter_Table * Table)
{
  int i;
  double x;

  double * P = Table->MEq->Probability_Distribution;
  int No_of_CONFIGURATIONAL_STATES  = Table->MEq->No_of_CONFIGURATIONAL_STATES;

  
  x = 0.0;
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) x += P[i];
  printf("Before Normalization, Total Probability is %g\n", x); 
  
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) P[i] /= x;

  x = 0.0;
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) x += P[i];
  printf("After Normalization, Total Probability is %g\n", x); 
}

