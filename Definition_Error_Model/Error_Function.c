#include <MODEL.h>

#define RANDOM (gsl_rng_uniform ( r ))

extern int MY_ERROR_HANDLER;

void
my_error_handler (const char * reason, const char * file, int line, int gsl_errno)
{

  gsl_stream_printf ("ERROR", file, line, reason);

  MY_ERROR_HANDLER = 1;

  fflush (stdout);

  fprintf (stderr, "My GSL-based error handler invoked.\n");

  fflush (stderr);

  //Print_Press_Key(1,0,".");
  //abort ();
}

double Measurement_Error_Function (gsl_rng * r, Parameter_Table * P,  double y)
{
	double x;

	gsl_function Density_Function;

	Density_Function.function = &my_Density_Function;
	Density_Function.params   = P;

	P->DETERMINISTIC_CASES = y;

	MY_ERROR_HANDLER = 0;

	gsl_error_handler_t * old_handler = gsl_set_error_handler ( &my_error_handler );

	x = da_gsl_ran_continuous_Function( r, &Density_Function, 0, 1);

	gsl_set_error_handler (old_handler);

	if (  MY_ERROR_HANDLER == 1 ) x = GSL_NAN;

	return( x );
}

double my_Density_Function (double x, void * params)
{
  double f, y;

  Parameter_Table * P = (Parameter_Table *) params;

  y = P->DETERMINISTIC_CASES;

  if( x < 0.0 && y < 0.0 )
    { f = 0.0; }
  else
    { f = 1.0 / (y * (2.0-exp(-1.0))) * exp( - fabs(x - y)/fabs (y) ) ; }

  return f;
}

double da_gsl_ran_continuous_Function( const gsl_rng * r, gsl_function * Density_Function,
				       int i, int N )
{
  /*
     This function returns a random integer from a
     continuous distribution with parameters as defined in structure
     Pam. Random numbers are generated by using, generically,
     the inversion method.

     For instance, the probability distribution for a power law
     distributed random variable, X, is:

     p(x) = P{ X = x } = C x^{-alpha}, where C is the normalization
     constant
  */
  double y, a;
  Parameter_Table * Pam = (Parameter_Table *)Density_Function->params;

  y = RANDOM;
  /* Application of the inversion method */
  /* Problem:

     Find a such that F(a) = y,

     where F(a) is the cummulative distribution function
     corresponding to the density function provided by
     gsl_function * my_Density_Function
  */
  const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
  static gsl_root_fsolver * s;

  if( i == 0 ){
    s            = gsl_root_fsolver_alloc (T);
#if defined VERBOSE
    printf ("using %s method\n", gsl_root_fsolver_name (s));
#endif
  }

  /* First Step: Determining the superior bracketing point x_upper */
  double I = 0.0;
  double x_upper = 100.0;
  double x_lower = 0.0;

  while( I <= y){

    x_upper = 2.0 * x_upper;

    I = Cummulative_Distribution_Function( x_upper, Density_Function );

  }

  /* Second step: Determining the point 'a' at give accuracy */
  int iter, max_iter, status;
  double x_lo, x_hi;
  Parameter_Table_Root_Solver P_R;
  gsl_function f_R;

  P_R.f = Density_Function;
  P_R.r = y;

  f_R.function = &Function_Root_Solver;
  f_R.params   = &P_R;

  gsl_root_fsolver_set (s, &f_R, x_lower, x_upper);

#if defined DA_DEBUGGING
  printf ("%5s [%9s, %9s] %9s %9s\n",
	  "iter", "lower", "upper", "root", "err(est)");
#endif

  iter = 0;
  max_iter = 100;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      a = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

#if defined DA_DEBUGGING
      if (status == GSL_SUCCESS)
	printf ("Converged:\n");

      printf ("%d [%.7f, %.7f] %.7f %.7f\n",
	      iter, x_lo, x_hi, a, x_hi - x_lo);
#endif
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  if( i == (N-1) ) gsl_root_fsolver_free ( s );

  return(a);
}

double Cummulative_Distribution_Function( double x, gsl_function * Density_Function )
{
  double I;
  double x_lo, x_hi;
  double Error;
  size_t Subinterval_Limit = 1.0e+4;
  int key                  = 4;
  double epsabs            = 1.0e-4;
  double epsrel            = 0.0;

  /*
     GSL_INTEG_GAUSS15  (key = 1)
     GSL_INTEG_GAUSS21  (key = 2)
     GSL_INTEG_GAUSS31  (key = 3)
     GSL_INTEG_GAUSS41  (key = 4)
     GSL_INTEG_GAUSS51  (key = 5)
     GSL_INTEG_GAUSS61  (key = 6)
   */

  Parameter_Table * P = (Parameter_Table *)Density_Function->params;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (Subinterval_Limit);

  x_lo = 0;
  x_hi = x;

  /* int status = da__gsl_integration_qag( Density_Function, x_lo, x_hi, */
  /* 					epsabs, epsrel, Subinterval_Limit, key, workspace, &I, &Error ); */

  int status = gsl_integration_qag( Density_Function, x_lo, x_hi,
				    epsabs, epsrel, Subinterval_Limit, key, workspace, &I, &Error );

  gsl_integration_workspace_free ( workspace );

  if( status != GSL_SUCCESS) I = GSL_NAN;

  return ( I );
}

double Function_Root_Solver( double x, void * p )
{
  double I;

  Parameter_Table_Root_Solver * Pam = (Parameter_Table_Root_Solver *)p;

  gsl_function * f = Pam->f;
  double         r = Pam->r;

  I = Cummulative_Distribution_Function( x, f ) - r;

  return( I );
}
