#include <MODEL.h>

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				      int * No_of_CONFIGURATIONAL_STATES,
				      int * n_DIMENSION,
				      int * n_x, int * n_y, int * n_z)
{
  Master_Equation * ME = Table->MEq; 
  
  * n_DIMENSION = 1;
  * n_z         = 1;
  * n_y         = 1;
  * n_x         = Table->TOTAL_No_of_CONSUMERS + 1;
  
  * No_of_CONFIGURATIONAL_STATES = Table->TOTAL_No_of_CONSUMERS + 1;
}

void Labels_for_Marginal_Probabilities (Parameter_Table * Table)
{
  Master_Equation * ME = Table->MEq; 
 
  /* Label of the Probability Dimensions */
  sprintf(ME->Marginal_Probability_Label[0], "Number of Free Consumers");

}

double Average_Number_of_Feeding_Events( Parameter_Table * Table, double T )
{
  double n_T;
  double n_R, m_0, K_R, Theta, Nu;

  K_R   = (double)Table->K_R;
  n_R   = (double)Table->TOTAL_No_of_RESOURCES; /* which determines resource density */
  m_0   = (double)Table->TOTAL_No_of_CONSUMERS;

  Theta = Table->Alpha_C_0 * n_R/K_R;  /* -H9  [Alpha_C_0] */
  Nu    = Table->Nu_C_0;               /* -H10 [Nu_C_0]    */
  
  n_T = Theta * Nu/(Nu+Theta) * ( T + Theta/Nu/( Nu + Theta) * (1.0-exp(-(Theta+Nu)*T)) ) * m_0;
  
  return (n_T);
}

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME )
{
  /* Since this process is one-dimensional, the probability distribution in vector
     form coincides with the same probability distribution... 
  */
  
  int i;
  int n,m,l;
  int s; 
  
  double * y = ME->Probability_Distribution; 

  assert(ME->n_DIMENSION == 1); 
  
  if (ME->n_DIMENSION == 1) 
    for(i=0; i<ME->n_x; i++) ME->P_n[i] = y[i];
  
}

void Marginal_Probability_Calculation ( Parameter_Table * Table )
{
  /* Since this process is one-dimensional, the marginal is the same as 
     probability distribution... 
  */ 
  int n, A_0;
  double S;
  
  Master_Equation * ME = Table->MEq;
  
  Probability_Distribution_Vector_into_Matrix_Form( ME );

  A_0 = Table->TOTAL_No_of_CONSUMERS;

  for( n = 0; n <= A_0; n++ )        /* Free Predators */
    ME->P_n_Marginal[n] = ME->Probability_Distribution[n];

}

void Marginal_Probability_Averages_Calculation ( Parameter_Table * Table )
{
  /* Since this process is one-dimensional, the marginal average is the same as 
     expected value of number of free predators based on the probability distribution. 
  */ 
  int n, m, A_0;
  double S;
  
  Master_Equation * ME = Table->MEq;
  
  Probability_Distribution_Vector_into_Matrix_Form( ME );

  A_0 = Table->TOTAL_No_of_CONSUMERS;

  S = 0.0;
  for( n = 0; n <= A_0; n++ )        /* Free Predators */  
      S += (double)n * ME->P_n[n];

  ME->Vector_Model_Variables[0] = S; 
}

void Print_Marginal_Averages( double Time_Current, Parameter_Table * Table)
{
  Master_Equation * ME = Table->MEq;
  
  assert(ME->n_DIMENSION == 1);
  
  printf("t = %g\t<n> = %g\n", Time_Current, ME->Vector_Model_Variables[0]);
}

void Print_Probability_Distribution ( Parameter_Table * Table )
{ 
  int n, m, A_0;
  
  Master_Equation * ME = Table->MEq;
  
  Probability_Distribution_Vector_into_Matrix_Form( ME );

  A_0 = Table->TOTAL_No_of_CONSUMERS;

  for( n = 0; n <= A_0; n++ )        /* Free Predators */
    printf("%.2g ", ME->P_n[n]);
    
  printf("\n"); 
  
}




