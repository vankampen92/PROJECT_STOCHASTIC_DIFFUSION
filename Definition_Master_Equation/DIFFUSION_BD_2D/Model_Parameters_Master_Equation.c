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
  // Print_Press_Key(1,0,"."); 
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

void Print_Marginal_Averages( double Time_Current, Parameter_Table * Table)
{
  Master_Equation * ME = Table->MEq;
  
  assert(ME->n_DIMENSION == 2);
  
  printf("t = %g\t<n> = %g\t<m> = %g\n",
	   Time_Current,
	   ME->Vector_Model_Variables[0],
	   ME->Vector_Model_Variables[1]);
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

void Saving_Marginal_Distribution_Triplets(Parameter_Table * Table, int j, double Time_Current)
{
  /* This function saves the triplets marginal probability distributions at a particular time, 
     Time_Current. 

     As an input argument, it takes Table. A member of Table is MEq, which stores all the 
     necessary information to save the whole distribution arising from the numerical integration
     of the master equation or related calculations, and, of course, the marginals.
     
     This function is specific of MODEL=DIFFUSION_BD_2D. That is the reason why it is stored
     in the corresponding model directory. 
   */
  /* Notice that the probability distribution is associated to a Time_Current around the 
     j-th time in Time_Vector[j] 
  */
  int i, l, n, m, A_0, ARA_MAX;
  double S;
  int No_of_POINTS;
  
  Master_Equation * ME = Table->MEq;
  
  double * y;
  double * x; 

  assert(Table->TYPE_of_MODEL == 13);             // MODEL=DIFFUSION_BD_2D

  A_0          = Table->TOTAL_No_of_CONSUMERS;    // A_0 + 1 = ME->n_x   
  assert(A_0%2 == 0);
  ARA_MAX      = Table->TOTAL_No_of_CONSUMERS/2;  
  No_of_POINTS = 1 + Table->TOTAL_No_of_CONSUMERS/2;  
  x = (double *)calloc(No_of_POINTS, sizeof(double));
  y = (double *)calloc(No_of_POINTS, sizeof(double));

  for(i=0; i < No_of_POINTS; i++) x[i] = (double)i;

  for (l=0; l < No_of_POINTS; l++) {
    S = 0.0;
    for(n=0; n < ME->n_x-2*l; n++)
      S += ME->P_nm[n][A_0-n-2*l];

    y[l] = S; 
  }
 
  char * Marginal = (char *)calloc(100, sizeof(char));
  char * pFile;

  pFile = strcat(Marginal, "Marginal_Probability_Triplets");
  pFile = strcat(Marginal, "_Time_");

  printf("Saving Marginal Probability for Triplets at Time %g\n", Time_Current);  
  Saving_to_File_double(Marginal, x, y, No_of_POINTS, j);  
  
  free(Marginal); 
  free(x); free(y); 
}




