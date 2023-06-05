#include <MODEL.h>

double PS_nn_Function( int n, int m, double C0, double p );
double PS_nk_Function( int n, int k, double Cn, double q ); 

void Stationary_Probability_Distribution (Parameter_Table * Table)
{
  double K_R, y_R, Theta;
  double Theory;
  double qI, q, p, C0,Cn;
  int    n_A, nA2;                     /* Total Number of Consumers */
  int n, k, m;

  Master_Equation * ME = Table->MEq;
  
  K_R   = (double)Table->K_R;
  y_R   = Table->TOTAL_No_of_RESOURCES; /* Resource density (in number of resource units) */
  Theta = Table->Alpha_C_0 * y_R/K_R;

  q  = Theta/Table->Nu_C_0;

  qI = Table->Nu_C_0/Theta;
  
  p  = Table->Eta_C_0/Table->Chi_C_0 * K_R;
 
  n_A = Table->TOTAL_No_of_CONSUMERS;

  assert(n_A%2 == 0);            /* Total Number of Consumers is an even number */   

  nA2 = n_A/2;                   /*                                             */
  
  C0 = ME->PS_nm[0][0] = 1.0;    /* Normalization Constant                      */
  
  for(n = 1; n<nA2; n++) {
    Cn = ME->PS_nm[n][n]  = PS_nn_Function(n, nA2, C0, p);
    
    for(k = 1; k<= n; k++) 
      ME->PS_nm[n+k][n-k] = PS_nk_Function(n, k, Cn, qI); 

    for(k = 1; k<= n; k++) 
      ME->PS_nm[n-k][n+k] = PS_nk_Function(n, k, Cn, q); 
  }

  /* Normalization                                                             */
  Norma_2D_Nx_Ny(ME->PS_nm, ME->n_x, ME->n_y);
  
  /* Calculation of the Theoretical Marginals with the normalized distribution */
  Marginal_Stationary_Probabilities_Calculation ( Table );
}

double PS_nn_Function( int n, int m, double C0, double p )
{
  int k; 
  double x;

  x = 0.0;
  for(k = 0; k < n; k++) 
    x += ( log((double)(m-k)) - 2.0*log((double)(n-k)) ); 

  x += (double)n * log(p) + log(C0); 
  
  x = exp(x); 
  
  return(x);
}

double PS_nk_Function( int n, int k, double Cn, double q )
{ 
  double x;
  int j;
  
  x = 0.0;
  for(j = 0; j < k; j++) 
    x += ( log((double)(n-j)) - log((double)(n+k-j)) ); 

  x += (double)k * log(q) + log(Cn); 
  
  x = exp(x); 
  
  return(x);
}

void Marginal_Stationary_Probabilities_Calculation ( Parameter_Table * Table )
{ 
  int l, n, m, A_0;
  int No_of_POINTS; 
  double S;
  
  Master_Equation * ME = Table->MEq;
  
  A_0 = Table->TOTAL_No_of_CONSUMERS;

  for( n = 0; n <= A_0; n++ ) {       /* Free Predators */
    S = 0.0;
    for( m = 0; m <= A_0 - n; m++ ) {
      S += ME->PS_nm[n][m];
    }
    ME->PS_n_Marginal[n] = S;
  }
  
  for( m = 0; m <= A_0; m++ ) {       /* Handling Predators */
    S = 0.0;
    for( n = 0; n <= A_0 - m; n++ ) {
      S += ME->PS_nm[n][m];
    }
    ME->PS_m_Marginal[m] = S;
  }
}





