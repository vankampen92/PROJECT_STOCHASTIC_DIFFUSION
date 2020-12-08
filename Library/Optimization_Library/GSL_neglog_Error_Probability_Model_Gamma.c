#include <MODEL.h>

double GSL_neglog_Error_Probability_Model_Gamma( double Data, double Theory,
						     Parameter_Fitting  * F )
{
  double neg_x, value_0, value_1;
  double p, a, b;
  int NOT_a_NUMBER; 
  int ILL_DEFINED_CALL; 

  Parameter_Table * Table = F->Table; 

  /* If X is a random variate following the gamma distribution, then: 
     
                             E[X]   = a*b;                          (1)
                             
                             Var[X] = a*b^2;                        (2)
  
     We assume here that the probability of a 'Data' value given a 'Theory' 
     prediction is a gamma distribution. We make two extra assumptions: 

     A. We assume that the expected value of a measure of 'Data' is well 
     predicted by our 'Theory' value. This involves that we believe we have
     a good average theory. Then:
     
     E[Data] = a * b = Theory;

     B. We assume that the standard error in 'Data' can be represented in 
     relative terms with respect the 'Theory' value, this is:
     
     Error[Data] = sqrt(Var[Data]) = epsilon; 

     These two conditions define parameters a and b in (1) and (2).
  */

  ILL_DEFINED_CALL = 0;  
  
  b = Table->Error_Parameter_0 * Table->Error_Parameter_0 / Theory;

  a = (Theory / Table->Error_Parameter_0) * (Theory / Table->Error_Parameter_0);
  
  if (Table->Error_Parameter_0 <= 0.0) ILL_DEFINED_CALL = 1;
  if (Theory                   <= 0.0) ILL_DEFINED_CALL = 1; 

  if (ILL_DEFINED_CALL == 0) {   
   
    // p = log( gsl_ran_gamma_pdf(Data, a, b) );

    // p = -gsl_sf_lngamma(a) -a*log(b) -Data/b + (a-1)*log(Data);

    assert(Data - Table->T->EPSILON > 0.0);

    value_0 = gsl_cdf_gamma_Q((Data + Table->T->EPSILON)/b, a, 1.0);
    value_1 = gsl_cdf_gamma_Q((Data - Table->T->EPSILON)/b, a, 1.0);

    /* value_0 = gammQ_double(a, (Data + Table->T->EPSILON)/b); */
    /* value_1 = gammQ_double(a, (Data - Table->T->EPSILON)/b); */

    p       = value_1 - value_0;
    
    NOT_a_NUMBER = 0; 
    NOT_a_NUMBER = gsl_isnan( p ); 
    
    if (NOT_a_NUMBER == 1) {
      printf(" p is fucking zero...\n");
      printf(" p = P{ x > Data-Epsilon } - P{ x > Data+Epsilon } = %g - %g\n", value_1, value_0);
      printf(" p is not a number in GSL_neglog_Error_Probability_Model(...)!!!\n");
      printf(" Theory = %g\t Data = %g\t Error = %g\n",
	     Theory, Data, Table->Error_Parameter_0);
      printf(" The Gamma Probability Density Function is Ill-defined:\n"); 
      printf(" NegLoglikelihood will be arbitrarily set to 0.0\n");
      printf(" Theory = %g\t Data = %g\t Error = %g\t a = %g\t b = %g\n",
	     Theory, Data, Table->Error_Parameter_0, a, b);
      exit(1);
    }
    else {
      if (p == 0.0) neg_x = 0.0;
      else {
	neg_x = -log(p);
	if (gsl_isnan( neg_x ) == 1) neg_x = 0.0;
      }
    }
    
  }
  else ILL_DEFINED_CALL = 1; 
  
  if( ILL_DEFINED_CALL ==  1) {
    neg_x = 0.0;
    printf(" Gamma Density is Ill-defined: Gamma(Data|Theory, Par) = %g\n NegLoglikelihood will be arbitrarily set to 0.0\n Theory = %g\t Data = %g\t Error = %g\t a = %g\t b = %g\n", p, Theory, Data, Table->Error_Parameter_0, a, b);
    
  }
    
  return(neg_x);
}
