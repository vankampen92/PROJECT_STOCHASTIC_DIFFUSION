#include <MODEL.h>

double GSL_neglog_Error_Probability_Model_Gaussian( double Data, double Theory,
						    Parameter_Fitting  * F )
{
  double neg_x, value_0, value_1;
  double p, x, sigma;
  int NOT_a_NUMBER;
  int ILL_DEFINED_CALL;

  Parameter_Table * Table = F->Table;
  /* If X is a Gaussian variate with zeron mian, then:

                             E[X]   = 0.0;                          (1)

                             Var[X] = sigma^2;                        (2)

     We assume here that the probability of the error in 'Data' around
     a given a 'Theory' value is Gaussian distributed with zero mean.
     prediction is a gamma distribution. We make two extra assumptions:

     We could also assume that the standard error in 'Data' can be represented in
     relative terms with respect the 'Theory' value, this is:

     Error[Data] = sqrt(Var[Data]) = epsilon * Theory;

     Therefore,
                            sigma = epsilon * Theory;

     Other assumptions are equally plausible.
  */
  F->Ill_Defined_Function_Call = 0; 
  ILL_DEFINED_CALL = 0;

  // if (Table->Error_Parameter_0 <= 0.0) ILL_DEFINED_CALL = 1;
  if (Table->Err_0 <= 0.0)             ILL_DEFINED_CALL = 1;
  if (Theory                   <= 0.0) ILL_DEFINED_CALL = 1;

  if (ILL_DEFINED_CALL == 0) {

    // sigma = Table->Error_Parameter_0 * Theory;
    // sigma   = Table->Error_Parameter_0;
    
    sigma   =  Table->Err_0; 
    
    x       = Data - Theory;
    
    value_0 = gsl_cdf_gaussian_Q(x + Table->T->EPSILON, sigma);
    value_1 = gsl_cdf_gaussian_Q(x - Table->T->EPSILON, sigma);

    p       = value_1 - value_0;

    NOT_a_NUMBER = 0;
    NOT_a_NUMBER = gsl_isnan( p );

    if (NOT_a_NUMBER == 1) {
      printf(" p is fucking zero or not even a fucking number!!!\n");
      printf(" p = P{ x > Data-Epsilon } - P{ x > Data+Epsilon } = %g - %g\n", value_1, value_0);
      printf(" p is not a number in GSL_neglog_Error_Probability_Model(...)!!!\n");
      printf(" Theory = %g\t Data = %g\t Error = %g\n",
	     Theory, Data, Table->Err_0);
      printf(" NegLoglikelihood will be arbitrarily set to a very large number \n");
      // exit(1);
      ILL_DEFINED_CALL = 1; 
    }

    else if (p <= 0.0 || p > 1) ILL_DEFINED_CALL = 1;
    
    else {
	neg_x = -log(p);
	if (gsl_isnan( neg_x ) == 1) ILL_DEFINED_CALL = 1;
    }
  }

  if( ILL_DEFINED_CALL ==  1) {
    neg_x = DBL_MAX;
    printf(" Gaussian Error Model is Ill-defined: NLL(Data|Theory, Par)\n NegLoglikelihood will be arbitrarily set to a very large\n Theory = %g\t Data = %g\t Error = %g\n", Theory, Data, Table->Err_0);

    F->Ill_Defined_Function_Call = 1; 
  }

  return(neg_x);
}
