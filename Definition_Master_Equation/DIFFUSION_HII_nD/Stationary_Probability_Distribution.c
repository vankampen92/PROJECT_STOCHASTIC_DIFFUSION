#include <MODEL.h>

void Multinomial_Stationarity_Probability_Distribution(double * p, int N, Parameter_Table * Table); 
void Binomial_Stationarity_Marginal_Probabilities(double p, int N, Parameter_Table * Table);

void Stationary_Probability_Distribution (Parameter_Table * Table)
{
  double Theory;
  double SR;
  int S, D, n, m, N, i;

  Master_Equation * ME = Table->MEq;
  double          * y = ME->Probability_Distribution; 
  
  S = Table->No_of_RESOURCES;
  D = ME->n_DIMENSION; 
  N = Table->TOTAL_No_of_CONSUMERS;

  assert(S == D);
  double * p = (double *)calloc(S+1, sizeof(double));

  SR = 0.0; 
  for (i = 0; i < S; i++)
    SR += (Table->Theta_Consumers[i]/Table->Nu_Consumers[i]);   
  
  for (i = 0; i < S; i++)
    p[i] = (Table->Theta_Consumers[i]/Table->Nu_Consumers[i])/(1.0 + SR);
  
  p[S] = 1.0/(1.0 + SR);
  
  Multinomial_Stationarity_Probability_Distribution(p, N, Table); 

  /* Normalization over the configurations space at D dimensions */
  D_Normalilzation_Probability_Distribution ( Table );

  if (ME->n_DIMENSION == 1) { 
    for( i = 0; i <= N; i++ ) 
      ME->PS_n[i] = y[i];
  }
  if (ME->n_DIMENSION == 2) { 
    for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++) {
      n = ME->Co[i]->n[0];
      m = ME->Co[i]->n[1];
      ME->PS_nm[n][m] = y[i];
    }
  }
  /* Calculation of the Theoretical Marginals with the normalized distribution */
  Marginal_Stationary_Probabilities_Calculation ( Table );
}

void D_Normalilzation_Probability_Distribution( Parameter_Table * Table )
{
  /* Normalization over the configuration space at D dimensions */  
  int i; 
  double S; 

  Master_Equation * ME = Table->MEq;
  double * y = ME->Probability_Distribution;                    
    
  S = 0.0;
  for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++)
    S += y[i]; 

  assert(S>0.0);

  for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++)
    y[i] /= S; 

  printf("Checking Normalization (full distribution). The sum S [%g] should be around 1\n", S); 
  Press_Key();
}

void Multinomial_Stationarity_Probability_Distribution(double * p, int N, Parameter_Table * Table)
{
  /* Input arguments: 
     . p is the vector of probabilities: (p_1, ..., p_S, p_{S+1})
     . N is the Table->TOTAL_No_of_CONSUMERS;

    Notice that this distribution becomes benomial when D = 1 (one resource type)

     Output: 
     . y = Table->MEq->Probability_Distributon[] is the multinomial probability at stationarity over
       configurations, from i=0 to No of CONFIGURATIONAL STATES-1.       
  */
  int S, D, i, k, n_H; 
  double X; 

  Master_Equation * ME = Table->MEq;  
  double * y = ME->Probability_Distribution; 
  int ** C   = ME->TabConfi; 

  S = Table->No_of_RESOURCES; 
  D = ME->n_DIMENSION;
  assert( S == D);

  int * n = (int *)calloc(D+1, sizeof(int));

  for(i=0; i < ME->No_of_CONFIGURATIONAL_STATES; i++) { 
    n_H = 0; 
    for(k = 0; k<D; k++) n_H += C[i][k]; 
     
    for(k = 0; k<D; k++) n[k] = C[i][k]; 
    n[S] = N - n_H; 
      X  = gsl_ran_multinomial_lnpdf(D+1, p, n);
    y[i] =  exp(X);     
  }

  free(n);    
}

void Marginal_Stationary_Probabilities_Calculation ( Parameter_Table * Table )
{ 
  /* The marginals can be calcualted as a binonmial theoretical distibution with parameters 
     N and p or by using the full multivariate multinomial distribution. Here, the full 
     distribution over configurational states is used.
  */
  int i, k, l, n, m, A_0;
  int No_of_POINTS; 
  double S;
  
  Master_Equation * ME = Table->MEq;
  A_0 = Table->TOTAL_No_of_CONSUMERS;
  double * y = ME->Probability_Distribution; 

  /* Maginal Probability Distributions for each (resource) Dimension k */
  for (k = 0; k<ME->n_DIMENSION; k++) {
    for (n = 0; n <= A_0; n++) {
      ME->MPD_S[k][n] = 0.0;
      for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++)
        if (ME->Co[i]->n[k] == n) 
              ME->MPD_S[k][n] += y[i];
    }
    
    /* Normalizing */
    S = 0.0; 
    for (n = 0; n<= A_0; n++) 
      S+= ME->MPD_S[k][n];

    for (n = 0; n<= A_0; n++) 
      ME->MPD_S[k][n] /= S;
      
    printf("Checking Normalization (marginals). The sum S [%g] should be around 1\n", S); 
    Press_Key();
  }   

  if (ME->n_DIMENSION == 1) {
    /* The marignal is the same as the full distribution */
    for( n = 0; n <= A_0; n++ ) {
      ME->PS_n_Marginal[n] = y[n];
      ME->PS_n[n]          = y[n];
    }
  }    
  if (ME->n_DIMENSION == 2) {
    for( n = 0; n <= A_0; n++ )        /* Handling Consumers on Resource 0 */
      ME->PS_n_Marginal[n] = ME->MPD_S[0][n];

    for( m = 0; m <= A_0; m++ )        /* Handling Consumers on Resource 1 */
      ME->PS_m_Marginal[m] = ME->MPD_S[1][m];
  }
}

void Binomial_Stationarity_Marginal_Probabilities(double p, int N, Parameter_Table * Table)
{
  /* The marginals can be calcualted as a binonmial theoretical distibution with parameters 
     N and p or by using the full multivariate distribution. Here, the first method is used
  */
  double SR; 
  int k,n, S, D;  
  
  Master_Equation * ME = Table->MEq;
  
  S = Table->No_of_RESOURCES; 
  D = ME->n_DIMENSION;
  assert( S == D);

  assert (N  == Table->TOTAL_No_of_CONSUMERS);
  SR = 0.0; 
  for (k = 0; k < S; k++)
    SR += (Table->Theta_Consumers[k]/Table->Nu_Consumers[k]);   

  for (k = 0; k<ME->n_DIMENSION; k++) {
    p = (Table->Theta_Consumers[k]/Table->Nu_Consumers[k])/(1.0 + SR);

    for (n = 0; n <= N; n++) 
      ME->MPD_S[k][n] = gsl_ran_binomial_pdf(n, p, N);   
  }
}



