#include <MODEL.h>

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				      int * No_of_CONFIGURATIONAL_STATES,
				      int * n_DIMENSION,
				      int * n_x, int * n_y, int * n_z)
{
  Master_Equation * ME = Table->MEq; 
  
  * n_DIMENSION = 2;
  * n_z         = 1;
  * n_y         = 1 + Table->TOTAL_No_of_CONSUMERS;  
  * n_x         = 1 + Table->TOTAL_No_of_CONSUMERS;
  
  * No_of_CONFIGURATIONAL_STATES = (Table->TOTAL_No_of_CONSUMERS + 1)*(Table->TOTAL_No_of_CONSUMERS + 2);
  * No_of_CONFIGURATIONAL_STATES = (* No_of_CONFIGURATIONAL_STATES) / 2;

  printf("No of CONFIGURATIONAL STATES = %d\n", * No_of_CONFIGURATIONAL_STATES);
  Press_Key(); 
}

void Labels_for_Marginal_Probabilities (Parameter_Table * Table)
{
  Master_Equation * ME = Table->MEq; 
 
  /* Label of the Probability Dimensions */
  sprintf(ME->Marginal_Probability_Label[0], "Number of Free Consumers");
  sprintf(ME->Marginal_Probability_Label[1], "Number of Handling Consumers");
}

/* Mapping of configuration states (n,m) into just one index (i) and vicevera */ 
void i_to_nm_Map(Parameter_Table * Table, int i, int * n_i, int * m_i)
{
  int m, a_0;
  double S;

  a_0 = Table->TOTAL_No_of_CONSUMERS;
  
  m = 0;
  S = Sumando(a_0, 1);

  while ( (i - S) >= 0.0 ) {
    m++;

    S = Sumando(a_0, m+1);
  }

  * n_i = i - Sumando(a_0, m);
  * m_i = m;
}

void nm_to_i_Map(Parameter_Table * Table, int * i, int n, int m)
{
  int a_0 = Table->TOTAL_No_of_CONSUMERS; 
  
  * i = n + Sumando(a_0, m); 
}

double Sumando(int a_0, int m)
{
  double R;

  R = (double)m/2.0 * (2.0*(a_0 + 1.0) - m + 1.0); 

  return(R); 
}

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME )
{

  /* y[]  ----> P_nm[][] */
  
  int i;
  int n,m;
  
  double * y = ME->Probability_Distribution;
  Parameter_Table * Table = (Parameter_Table *)ME->Table; 

  int No_of_CONFIGURATIONAL_STATES = ME->No_of_CONFIGURATIONAL_STATES;
   
  for(i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {
    
    i_to_nm_Map(Table, i, &n, &m);   

    ME->P_nm[n][m] = y[i];
  }
}

void Marginal_Probability_Calculation ( Parameter_Table * Table )
{ 
  int n, m, A_0;
  double S;
  
  Master_Equation * ME = Table->MEq;
  
  Probability_Distribution_Vector_into_Matrix_Form( ME );

  A_0 = Table->TOTAL_No_of_CONSUMERS;

  for( n = 0; n <= A_0; n++ ) {       /* Free Predators */
    S = 0.0;
    for( m = 0; m <= A_0 - n; m++ ) {
      S += ME->P_nm[n][m];
    }
    ME->P_n_Marginal[n] = S;
  }

  for( m = 0; m <= A_0; m++ ) {       /* Handling Predators */
    S = 0.0;
    for( n = 0; n <= A_0 - m; n++ ) {
      S += ME->P_nm[n][m];
    }
    ME->P_m_Marginal[m] = S;
  }
}

void Marginal_Probability_Averages_Calculation ( Parameter_Table * Table )
{ 
  int n, m, A_0;
  double S;
  
  Master_Equation * ME = Table->MEq;
  
  Probability_Distribution_Vector_into_Matrix_Form( ME );

  A_0 = Table->TOTAL_No_of_CONSUMERS;

  S = 0.0;
  for( n = 0; n <= A_0; n++ ) {       /* Free Predators */
    for( m = 0; m <= A_0 - n; m++ ) {
      S += (double)n * ME->P_nm[n][m];
    }
  }

  ME->Vector_Model_Variables[0] = S; 

  S = 0.0;
  for( m = 0; m <= A_0; m++ ) {       /* Handling Predators */
    for( n = 0; n <= A_0 - m; n++ ) {
      S += (double)m * ME->P_nm[n][m];
    }
  }

  ME->Vector_Model_Variables[1] = S; 
}

void Print_Probability_Distribution ( Parameter_Table * Table )
{ 
  int n, m, A_0;
  
  Master_Equation * ME = Table->MEq;
  
  Probability_Distribution_Vector_into_Matrix_Form( ME );

  A_0 = Table->TOTAL_No_of_CONSUMERS;

  for( n = 0; n <= A_0; n++ ) {       /* Free Predators */
    for( m = 0; m <= A_0 - n; m++ ) {
      printf("%.2g ", ME->P_nm[n][m]);
    }
    printf("\n"); 
  }
 
}





