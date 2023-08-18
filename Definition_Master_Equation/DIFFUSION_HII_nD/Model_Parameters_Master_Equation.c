#include <MODEL.h>

void assert_HLL_nD_Equal_Nus( Parameter_Table * Table )
{
  int k; 
  int i; 
  int n_DIMENSION = Table->MEq->n_DIMENSION; 

  #if defined DIFFUSION_HII_nD
    printf(" The mathematical expression for the evolving probability for each configurational state\n");
    printf(" is only valid is handling times are the same for all resources\n");
    
    k = 0; 
    for(i=0; i<n_DIMENSION; i++)
      if (Table->Nu_Consumers[i] != Table->Nu_C_0) 
        k = 1; 

    if ( k == 1) {
      printf(" Some Handing Times differ.\n");
      printf(" Theoretical Evolving probabilities are not exact. Therefore, at most, they represent an approximation\n");
      printf("\n");
    }

    assert(k == 0);
  #endif
}

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				                              int * No_of_CONFIGURATIONAL_STATES,
				                              int * n_DIMENSION,
				                              int * n_x, int * n_y, int * n_z)
{
  /* Here, I maintain n_y and n_z, for compatibility, but I only use and set up n_x 
     as the vector defining the size of each dimension of the configurational state. 
     This vector should have been allocated before calling this function. 
  */
  int i; 
  Master_Equation * ME = Table->MEq; 
  
  * n_DIMENSION = Table->No_of_RESOURCES;
  * n_x         = 1 + Table->TOTAL_No_of_CONSUMERS;
  * n_z         = 1 + Table->TOTAL_No_of_CONSUMERS;
  * n_y         = 1 + Table->TOTAL_No_of_CONSUMERS;
  
  * No_of_CONFIGURATIONAL_STATES = Func_Hyper( Table->No_of_RESOURCES, Table->TOTAL_No_of_CONSUMERS );

  printf(" No of CONFIGURATIONAL STATES = %d\n", * No_of_CONFIGURATIONAL_STATES);
  printf(" TOTAL No of CONSUMERS        = %d\n", * n_DIMENSION);
  // Print_Press_Key(1,0,"."); 
}

void Labels_for_Marginal_Probabilities (Parameter_Table * Table)
{
  int i; 
  Master_Equation * ME = Table->MEq; 
 
  /* Label of the Probability Dimensions */
  for(i=0; i<Table->No_of_RESOURCES; i++)
    sprintf(ME->Marginal_Probability_Label[i], "Handling Consumers on the %d-Resource", i);
  
}

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME )
{

  /* y[]  ----> P[]...[] */
  
  int i;
  int * n = (int *)calloc(ME->n_DIMENSION, sizeof(int));

  configuration ** Co = ME->Co;
  int ** Configuration_Table = ME->TabConfi; 

  double * y = ME->Probability_Distribution;
  // Parameter_Table * Table = (Parameter_Table *)ME->Table; 

  int No_of_CONFIGURATIONAL_STATES = ME->No_of_CONFIGURATIONAL_STATES;
   
  for(i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {
    
    i_to_Configuration_Map(Configuration_Table, i, n, ME->n_DIMENSION);

    assert( ME->n_DIMENSION > 1 && ME->n_DIMENSION <= ME_n_DIMENSION_MAXIMUM);

    if (ME->n_DIMENSION == 2)  ME->P_nm[n[0]][n[1]] = y[i];
    if (ME->n_DIMENSION == 3)  ME->P_nml[n[0]][n[1]][n[2]] = y[i];
    // if (ME->n_DIMENSION == 4)  ME->P_nmlk[n[0]][n[1]][n[2]][n[3]] = y[i];
    // if (ME->n_DIMENSION == 5)  ME->P_nmlkj[n[0]][n[1]][n[2]][n[3]][n[4]] = y[i];
    // if (ME->n_DIMENSION == 6)  ME->P_nmlkji[n[0]][n[1]][n[2]][n[3]][n[4]][n[5]] = y[i];
    // if (ME->n_DIMENSION == 7)  ME->P_nmlkjih[n[0]][n[1]][n[2]][n[3]][n[4]][n[5]][n[6]] = y[i];
    // if (ME->n_DIMENSION == 8)  ME->P_nmlkjihg[n[0]][n[1]][n[2]][n[3]][n[4]][n[5]][n[6]][n[7]] = y[i];    
    // if (ME->n_DIMENSION == 9)  ME->P_nmlkjihgf[n[0]][n[1]][n[2]][n[3]][n[4]][n[5]][n[6]][n[7]][n[8]] = y[i];
    // if (ME->n_DIMENSION == 10) ME->P_nmlkjihgfe[n[0]][n[1]][n[2]][n[3]][n[4]][n[5]][n[6]][n[7]][n[8]][n[9]] = y[i];
    // Que locura de amor...
    else {
      printf(" Dimension of the Probability Distribution is bigger than 3!!! ");
      printf(" n_DIMENSION = %d", ME->n_DIMENSION);
      printf(" The probability distribution is not allocated as n_DIMENSION object (tensor)\n");
      printf(" Something went wrong... The program should not reach this point\n");
      printf(" The program will exit\n");
      exit(0);
    }
  } 

  free(n);
}

void Marginal_Probability_Calculation ( Parameter_Table * Table )
{ 
  int i, k, n, m, A_0;
  double S;
  
  Master_Equation * ME = Table->MEq;  
  A_0 = Table->TOTAL_No_of_CONSUMERS;

  if( ME->n_DIMENSION == 1) {
    for( n = 0; n <= A_0; n++ ) { 
      ME->P_n_Marginal[n] = ME->P_n[n];
    }
  }
  else if( ME->n_DIMENSION == 2) {
    Probability_Distribution_Vector_into_Matrix_Form( ME );

    for( n = 0; n <= A_0; n++ ) {       /* Fist Dimension: No of Handling Consumers */
      S = 0.0;                          /* on first resource [0] */
      for( m = 0; m <= A_0 - n; m++ ) {
        S += ME->P_nm[n][m];
      }
      ME->P_n_Marginal[n] = S;
    }

    for( m = 0; m <= A_0; m++ ) {       /* Second Dimension:  No of Handling Consumers */
      S = 0.0;                          /* on second resource [1] */
      for( n = 0; n <= A_0 - m; n++ ) {
        S += ME->P_nm[n][m];
      }
      ME->P_m_Marginal[m] = S;
    }
  }
  else {
    double * y = ME->Probability_Distribution; 

    /* Maginal Probability Distributions for each (resource) Dimension k */
    for (k = 0; k<ME->n_DIMENSION; k++) {
      for (n = 0; n<= A_0; n++) {
        ME->MPD[k][n] = 0.0;
        for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++)
          if (ME->Co[i]->n[k] == n) 
            ME->MPD[k][n] += y[i];
      }
      /* Normalizing */
      S = 0.0; 
      for (n = 0; n<= A_0; n++) 
        S+= ME->MPD[k][n];
      for (n = 0; n<= A_0; n++) 
        ME->MPD[k][n] /= S;
    }   
  } 
}

void Marginal_Probability_Averages_Calculation ( Parameter_Table * Table )
{ 
  int i, n, m, k, A_0;
  double S;
  
  Master_Equation * ME = Table->MEq;
  A_0 = Table->TOTAL_No_of_CONSUMERS;

  assert (Table->No_of_RESOURCES == ME->n_DIMENSION );
  
  if( ME->n_DIMENSION == 2) {
    Probability_Distribution_Vector_into_Matrix_Form( ME );
  
    S = 0.0;
    for( n = 0; n <= A_0; n++ ) {        /* Fist Dimension: No of Handling Consumers */
      for( m = 0; m <= A_0 - n; m++ ) {  /* on first resource [0] */
        S += (double)n * ME->P_nm[n][m];
     }
    }
    ME->Vector_Model_Variables[0] = S; 

    S = 0.0;
    for( m = 0; m <= A_0; m++ ) {       /* Second Dimension:  No of Handling Consumers */
      for( n = 0; n <= A_0 - m; n++ ) { /* on second resource [1] */
        S += (double)m * ME->P_nm[n][m];
      }
    }
    ME->Vector_Model_Variables[1] = S;
  }
  else {
    double * y = ME->Probability_Distribution; 

    /* Maginal Probability Distributions for each Dimension k */
    for (k = 0; k<ME->n_DIMENSION; k++) {
      S = 0.0; 
      for (n = 0; n<= A_0; n++) {
        for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++) {
          if (ME->Co[i]->n[k] == n) S += (double)n * y[i];
        }
      }
      ME->Vector_Model_Variables[k] = S;
    }
  }
}

void Print_Marginal_Averages( double Time_Current, Parameter_Table * Table)
{
  int k; 
  Master_Equation * ME = Table->MEq;
  
  if(ME->n_DIMENSION == 2) {
  
    printf("t = %g\t<n> = %g\t<m> = %g\n",
	          Time_Current,
	          ME->Vector_Model_Variables[0],
	          ME->Vector_Model_Variables[1]);
  }
  else {

    printf("t = %g\t", Time_Current);
    for (k = 0; k<ME->n_DIMENSION; k++) 
      printf("<n>[%d] = %g\t", k, ME->Vector_Model_Variables[k]);
    
    printf("\n");
  } 

  Print_Press_Key(1,0,"."); 
}

void Print_Probability_Distribution ( Parameter_Table * Table )
{ 
  int i, n, m, A_0;
  
  Master_Equation * ME = Table->MEq;
  A_0 = Table->TOTAL_No_of_CONSUMERS;
 
  if (ME->n_DIMENSION <= 2) {

    Probability_Distribution_Vector_into_Matrix_Form( ME );

    for( n = 0; n <= A_0; n++ ) {
      for( m = 0; m <= A_0 - n; m++ ) {
        printf("%.2g ", ME->P_nm[n][m]);
      }
      printf("\n"); 
    } 
  }
  else {
    double * y = ME->Probability_Distribution;
    
    for(i=0; i<ME->No_of_CONFIGURATIONAL_STATES; i++) {
      printf(" P( c(%d)=[ ", i);
      for(m=0; m<ME->n_DIMENSION; m++) {
        printf("%d  ", ME->Co[i]->n[m]);
      }
      printf("] ) = %g\n", y[i]);
    }
  }
}

void Saving_Marginal_Distributions(Parameter_Table * Table, int j, double Time_Current)
{
  /* This function saves the marginal probability distributions at a particular time, 
     Time_Current. 

     As an input argument, it takes Table. A member of Table is MEq, which stores all the 
     necessary information to save the whole distribution arising from the numerical integration
     of the master equation or related calculations, and, of course, the marginals.
     
     This function is specific of MODEL=DIFFUSION_HII_nD. That is the reason why it is stored
     in the corresponding model directory. 
   */
  /* Notice that the probability distribution is associated to a Time_Current around the 
     j-th time in Time_Vector[j] 
  */
  int i, k, n, m, A_0;
  double S;
  int No_of_POINTS;
  char * Marginal  = (char *)calloc(100, sizeof(char));
  char * Dimension = (char *)calloc(10, sizeof(char));
  char * pFile;

  Master_Equation * ME = Table->MEq;
  A_0          = Table->TOTAL_No_of_CONSUMERS;    
  
  No_of_POINTS = 1 + A_0;  
  double * x = (double *)calloc(No_of_POINTS, sizeof(double));
  double * y = (double *)calloc(No_of_POINTS, sizeof(double));

  for(i=0; i < No_of_POINTS; i++) x[i] = (double)i;

  if(ME->n_DIMENSION == 2) {
  
    for (n=0; n < No_of_POINTS; n++) 
      y[n] = ME->P_n_Marginal[n]; 
    
    Marginal[0] = '\0';
    pFile = strcat(Marginal, "Marginal_Probability_0");
    pFile = strcat(Marginal, "_Time_");
    Saving_to_File_double(Marginal, x, y, No_of_POINTS, j);
    printf("%s: Marginal Probability (Dimension = 0) at Time %g saved!!!\n", Marginal, Time_Current);  
    
    for (n=0; n < No_of_POINTS; n++) 
      y[n] = ME->P_m_Marginal[n]; 
    
    Marginal[0] = '\0';
    pFile = strcat(Marginal, "Marginal_Probability_1");
    pFile = strcat(Marginal, "_Time_");  
    Saving_to_File_double(Marginal, x, y, No_of_POINTS, j);
    printf("%s: Marginal Probability (Dimension = 1) at Time %g saved!!!\n", Marginal, Time_Current);  
  }
  else {
    assert(ME->n_DIMENSION == Table->No_of_RESOURCES);

    for(k = 0; k<ME->n_DIMENSION; k++) { 
       for (n=0; n < No_of_POINTS; n++) 
        y[n] = ME->MPD[k][n];  
    
      Marginal[0] = '\0';
      pFile = strcat(Marginal, "Marginal_Probability_");
      sprintf(Dimension, "%d", k);
      pFile = strcat(Marginal, Dimension);  
      pFile = strcat(Marginal, "_Time_");  
      Saving_to_File_double(Marginal, x, y, No_of_POINTS, j);
      printf("%s: Marginal Probability (Dimension = %d) at Time %g saved!!!\n", Marginal, k, Time_Current);
    }
  }

  free(Marginal);
  free(Dimension); 
  free(x); free(y); 
}

/* Auxiliary functions (to the end of the file) to define and deal with 
   the space of configurations for DIFFUSION_HII_nD model. 

   This functions are fully model-specific!!! 

   This is why they lie on the model-specific directory 
*/
unsigned long long int Func_Hyper(unsigned int D, unsigned int k)
{
  /* Exhaustive count of the number of configuations.
     
     For instance, for D=3, configurations, (n_x, n_y, n_z), should fulfill 
     the condition: 
                             n_x + n_y + n_z <= k,
     where n_x, n_y, and n_z can take values between 0 and k.
     
     This function takes adventage of the recurrent relationship: 

                             S(D, k) = \sum_{j=0}^k S(D-1, j)
  */
  unsigned long long int S, i; 
  
  if(D == 0) return 1; 

  S = 0; 
  for(i=0; i<=k; i++)
    S += Func_Hyper(D-1, i);  
  
  return S;
}

void Groups(unsigned int D, int D_Max, int * V, int * s)
{
  /* This function calculates group configurations. 

     Vector V[] should be declared with enough room to host all groups 
  */
  int i, j, k; 
  int n_0, n; 

  if (D < D_Max ) {
    n_0 = V[0]; 

    j=0;
    k=0; 
    while ( n_0 > 0 )
    {
      for(i=0; i < *s; i++) {
        n = V[i]-k;
        if(n > 0)  
          V[j++] = n;  
      }
      n_0--;
      k++; 
    }
    *s = j;

    Groups(D+1, D_Max, V, s);
  }
}

unsigned long long int General_Configuration_Table(unsigned int D, unsigned int n_0, 
                                                   int ** Configuration_Table)
{
  /* A recursion to build exhaustively all possible configurations such that, for instance, 
     for D=3, configurations, (n_x, n_y, n_z), fulfill the condition:

                             n_x + n_y + n_z <= n_0,

     where n_x, n_y, and n_z can take values between 0 and n_0

    Input: 
      . D,    dimension. 
      . n_0,  condition: n_x + n_y + n_z <= n_0, (for instance, for D=3), 
      . Configuration_Table[][] has been allocated before calling this function

    Output: 
      . Configuraiton_Table[][] is filled with all possible configurations. 
  */
  int i, j, k, l, s, s_0, d, g_D;
  unsigned long long int S;
  int * G_c; 

  int ***Configuration  = (int ***)calloc(D, sizeof(int **));

  int * S_D             = (int *)calloc(D, sizeof(int));
  for(i=0; i<D; i++) 
    S_D[i] = Func_Hyper(i+1, n_0); // S_1(n_0), ..., S_D(n_0) //
  
  for(i=0; i<D; i++){
    Configuration[i] = (int **)calloc(S_D[i], sizeof(int *)); 
    for(j=0; j < S_D[i]; j++)
      Configuration[i][j] = (int *)calloc(i+1, sizeof(int));
  }

  /* Index i=0 in Configuration[i][][] corresponds to the initial configuration, d=1.   */
  /* Index i=1 in Configuration[i][][] corresponds to the configuration, d=2, and so on */
  /* Configuration for d=1 is set up here: */
  for(i=0; i<=n_0; i++) 
    Configuration[0][i][0] = i;

  /* After each round of the loop the next configuration, d+1, is calculated 
     based on the previous one */
  d=1; 
  while (d < D)
  {
    g_D = Func_Hyper(d-1, n_0);
    G_c = (int *)calloc(g_D, sizeof(int));
    G_c[0] = n_0+1; 
    s = 1; 
    Groups(1, d, G_c, &s);

    S=0; 
    for(i=0; i<=n_0; i++) {
      s_0=0;
      for(k=0; k<s ; k++) {
        for(j=0; j<G_c[k]-i; j++) {
            for(l=0; l < d; l++){ 
              Configuration[d][S][l] = Configuration[d-1][j+s_0][l];
            }
            Configuration[d][S][d] = i; 
            S++;
        }
        s_0 += G_c[k];  
      }
    }  

    free(G_c);
    assert(S_D[d] == S);
    d++;
  }

  /* Output: Storing final Configuration Table at D dimensions: */
  for(i=0; i < S; i++) 
     for(j=0; j < D; j++)
       Configuration_Table[i][j] = Configuration[D-1][i][j];

  for(i=0; i<D; i++){ 
    for(j=0; j < S_D[i]; j++)
      free(Configuration[i][j]);
    free(Configuration[i]);
  }
  free(Configuration);

  free(S_D);
  return(S);
}

/* Mapping of configuration states (n_1,..., n_D) into just one index (i) */ 
int Configuration_to_i_Map(int * Configuration, int D, int n_0)
{
  /* Mapping of configuration states (n_1,..., n_D) into just one index (i) */ 
  int i, k, ks, n, m; 

  i = 0; ks=0;
  for (n=D-1; n >= 0; n--) {
    k  = Configuration[n];

    for (m = 0; m < k; m++)
      i += Func_Hyper(n, n_0-ks-m);
    
    ks+= k;
  }

  return i;
}

void i_to_Configuration_Map(int ** Configuration_Table, int i, int * Configuration, int D)
{
  /* Selecting the i-th configuration state from Configuration Table and storing it in 
     Configuration array (n_1,..., n_D) 
  */ 
  int j;
  for(j=0; j < D; j++)
       Configuration[j] = Configuration_Table[i][j];
}

void Writing_General_Configuration_Table(int D,
                                         unsigned long long int S, 
                                         int ** Configuration_Table)
{
  /* S: Number of configuations... */
  int i, j;
  
  for(i=0; i < S; i++) {
    printf("Configuration[%d] = [", i);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table[i][j]);
    printf("]\n");
    if (i < 20 || i > S-20) getchar();
  }      
}      

void Comparing_Two_Configuration_Tables(int D, unsigned long long int S, 
                                        int ** Configuration_Table_D5, 
                                        int ** Configuration_Table)
{
  int i, j; 

  printf("Comparting two Configuration Tables...\n");
  
  for(i=0; i < S; i++) {
    printf("Configuration[%d] = [", i);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table_D5[i][j]);
    printf("]\t");

    printf("Configuration[%d] = [", i);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table[i][j]);
    printf("]\n");
    if (i < 20 || i > S-20) getchar();
  }
}

void Pickup_a_Configuration_from_Table(int D,
                                       int Index, 
                                       int ** Configuration_Table)
{
  /* Index: Order for a configuration... */
  int j;
  
  printf("Configuration[%d] = [", Index);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table[Index][j]);
  printf("]\n");
  // getchar();
}        

int Vector_of_Adjacent_Configurations(int * Configuration, int D, int n_0, int d,  
                                      int ** Nei, int * Index_Nei)
{
  int i, j, n, a; 

  n=0; 
  for(i=0; i<D; i++) 
    n += Configuration[i];   /* Controling Hyperplane Condtion!!! */
  
  a = 0;
  for(i=0; i<D; i++) {
     
    if( (Configuration[i]+d) >= 0 && (n+d) <= n_0 ) {
      for(j=0; j<D; j++) {
        if (j == i) Nei[a][j] = Configuration[j] + d;
        else        Nei[a][j] = Configuration[j];      
      }  
      a++;
    }  
  }

  for(i=0; i<a; i++) {
    Index_Nei[i] = Configuration_to_i_Map(Nei[i], D, n_0);
  }
  
  /* No of Neighboring configurations */
  return(a);
}

void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0)
  {
    int i, j, a, m, Index; 

    int ** Nei = (int **)calloc(D, sizeof(int *));
    for(i=0; i<D; i++)
      Nei[i] = (int *)calloc(D, sizeof(int));

    for (i=0; i<S; i++) {
      
          // printf(" P( c(%d)=[ ", i);
          for(m=0; m < D; m++) {
            Co[i]->n[m] = Configuration_Table[i][m];
            // printf("%d  ", Co[i]->n[m]);
          }
          // printf(" ] ) = 0.0\n"); 

      a = Vector_of_Adjacent_Configurations(Co[i]->n, D, n_0, +1, Nei, Co[i]->anUp);
      Co[i]->anUp[D] = a;

#if defined VERBOSE
      printf(" There are %d neighboring configurations around c[%d]=[ ", a, i);
      for(j=0; j<D; j++) printf("%d ", Co[i]->n[j]);   /* Configuration[i] = [n_1, ..., n_S] */
      printf("] transitioning down into it:\n");
#endif       
      for(j=0; j<a; j++){
        Index = Co[i]->anUp[j];
        // Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
      }

      a = Vector_of_Adjacent_Configurations(Co[i]->n, D, n_0, -1, Nei, Co[i]->anDw);
      Co[i]->anDw[D] = a; 

#if defined VERBOSE
      printf(" There are %d neighboring configurations around c[%d]=[ ", a, i);
      for(j=0; j<D; j++) printf("%d ", Co[i]->n[j]);   /* Configuration[i] = [n_1, ..., n_S] */
      printf("] transitioning up into it:\n");
#endif
      for(j=0; j<a; j++){
        Index = Co[i]->anDw[j];
        // Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
      } 

      // printf("...\n\n");
    }

    for(i=0; i<D; i++)
      free(Nei[i]);
    free(Nei);
  }

