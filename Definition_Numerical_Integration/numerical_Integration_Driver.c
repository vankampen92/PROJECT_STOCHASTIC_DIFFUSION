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

int numerical_Integration_Driver( Parameter_Table * Table, 
                                  int j, double * Time_Current )
{
  int status;

#ifdef GSL_BSIMP
  const gsl_odeiv_step_type * T = gsl_odeiv_step_bsimp; /* This 'bsimp' uses the Jacobian */
#elif defined GSL_RK4
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
#else
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; /* Default Numerical Integrator   */
#endif
  
  int MODEL_STATE_VARIABLES = Table->MODEL_STATE_VARIABLES;
  //Parameter_Model * P = Table->P;

  gsl_odeiv_step * s 
         = gsl_odeiv_step_alloc (T, MODEL_STATE_VARIABLES);
  gsl_odeiv_control * c 
         = gsl_odeiv_control_y_new (GSL_ODE_ERR_0, GSL_ODE_ERR_1);
  gsl_odeiv_evolve * e 
         = gsl_odeiv_evolve_alloc (MODEL_STATE_VARIABLES);
     
  gsl_odeiv_system sys = {function, jacobian, MODEL_STATE_VARIABLES, Table};
     
  double * t = Time_Current;
  double t1  = Table->T->Time_Vector[j];

  double h = 1e-6;
  
  double * y  = Table->Vector_Model_Variables;
    
  while (   (*t) < t1   )
    {
      status = gsl_odeiv_evolve_apply (e, c, s,
                                       &sys, 
                                       t, t1,
                                       &h, y);
     
      if (status != GSL_SUCCESS) break;

      #if defined VERBOSE
          if (Table->MODEL_STATE_VARIABLES >= 4 )
	          printf ("t = %.5g (%.5g, %.5g, %.5g, %.5g)\n", *t, y[0], y[1], y[2], y[3]);
      #endif
    }
     
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  
  return status;
}




