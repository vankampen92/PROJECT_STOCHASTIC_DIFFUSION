#include <MODEL.h>

void Stationary_Probability_Distribution (Parameter_Table * Table)
{
  
  double K_R, y_R, Theta;
  double Theory;
  double p;
  int    n; 
  unsigned int k;
  
  K_R   = (double)Table->K_R;
  y_R   = Table->TOTAL_No_of_RESOURCES; /* Resource density (in number of resource units) */
  Theta = Table->Alpha_C_0 * y_R/K_R;
  
  p = Table->Nu_C_0 / (Table->Nu_C_0 + Theta) ;

  n = Table->TOTAL_No_of_CONSUMERS;
  
  for(k = 0; k<= n; k++) {
    Table->MEq->PS_n[k] = gsl_ran_binomial_pdf(k, p, n);
    Table->MEq->PS_n_Marginal[k] = gsl_ran_binomial_pdf(k, p, n);
  }
}

