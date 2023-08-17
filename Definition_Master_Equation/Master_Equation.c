#include <MODEL.h>
/* 
   This functions allocate, initialize and free the structure 
   to conduct de numerical integration of a Master Equation 
*/
void Master_Equation_Allocation ( Master_Equation * ME,
				                          int No_of_CONFIGURATIONAL_STATES,
				                          int n_DIMENSION,
				                          int n_x,
				                          int n_y,
				                          int n_z )
{
  /* Input arguments: 
     
     . n_DIMENSION, dimension of the mutivariate probability distribution. This is can 
     be maximum 3, this is, P(n, m, l) 
   
     . No_of_CONFIGURATIONAL_STATES, which is equal to (or less than ) the product of 
     n_x * n_y * n_z
     
     . n_x, n_y, n_z, are the maximum sizes for each dimension. 

     If the natural support of P(n, m, l) does not go from n=0 to n_x, m=0 to n_y, and 
     l=0 to n_z, then P(n, m, l) will take zero values out of its actual support (see 
     calloc() calls for zero allocations below). For the sake of saving memmory, when 
     n_DIMENSION is bigger than 3, the probability distribution are no longer stored as 
     tensor objects.    
  */ 
  int i, j, k;
  Parameter_Table * Table = (Parameter_Table *)ME->Table; 

  int n_Time = Table->T->I_Time; 
  
  assert(No_of_CONFIGURATIONAL_STATES < MAX_No_of_CONFIGURATIONAL_STATES);
  
  ME->Evolving_Distribution = (double **)calloc(n_Time, sizeof(double *));
  for(k=0; k<n_Time; k++)
    ME->Evolving_Distribution[k]=(double *)calloc(No_of_CONFIGURATIONAL_STATES, sizeof(double));
  
  ME->Probability_Distribution=(double *)calloc(No_of_CONFIGURATIONAL_STATES, sizeof(double));

  ME->Probability_Distribution_Time_0 = (double *)calloc(No_of_CONFIGURATIONAL_STATES,sizeof(double));

  ME->MPD = (double **)calloc(n_DIMENSION, sizeof(double *));
    for(i=0; i<n_DIMENSION; i++)
      ME->MPD[i] = (double *)calloc(Table->TOTAL_No_of_CONSUMERS+1, sizeof(double));

  ME->MPD_T = (double ***)calloc(n_Time, sizeof(double **));
  for(k=0; k<n_Time; k++) {
    ME->MPD_T[k] = (double **)calloc(n_DIMENSION, sizeof(double *));
    for(i=0; i<n_DIMENSION; i++)
      ME->MPD_T[k][i] = (double *)calloc(Table->TOTAL_No_of_CONSUMERS+1, sizeof(double));
  }

  ME->MPD_S = (double **)calloc(n_DIMENSION, sizeof(double *));
    for(i=0; i<n_DIMENSION; i++)
      ME->MPD_S[i] = (double *)calloc(Table->TOTAL_No_of_CONSUMERS+1, sizeof(double));

  if (n_DIMENSION == 1) {
    ME->PS_n = (double *)calloc(n_x, sizeof(double) );       /* Stationary Distribution */
    ME->P_n  = (double *)calloc(n_x, sizeof(double) );
    ME->P_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->PS_n_Marginal = (double *)calloc(n_x, sizeof(double));
  }
  else if (n_DIMENSION == 2) {
    ME->PS_nm = (double **)calloc(n_x, sizeof(double *));    /* Stationary Distribution */
    ME->P_nm = (double **)calloc(n_x, sizeof(double *)); 
    for(i=0; i<n_x; i++) {
      ME->P_nm[i]  = (double *)calloc(n_y, sizeof(double));
      ME->PS_nm[i] = (double *)calloc(n_x, sizeof(double));  /* Stationary Distribution */
    }
    ME->P_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->P_m_Marginal = (double *)calloc(n_y, sizeof(double));

    ME->PS_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->PS_m_Marginal = (double *)calloc(n_y, sizeof(double));
  }
  else if (n_DIMENSION == 3) {
    ME->P_nml = (double ***)calloc(n_x, sizeof(double **)); 
    ME->PS_nml = (double ***)calloc(n_x, sizeof(double **));        /* Stationary Distribution */
    for(i=0; i<n_x; i++) {
      ME->P_nml[i]  = (double **)calloc(n_y, sizeof(double *));
      ME->PS_nml[i] = (double **)calloc(n_y, sizeof(double *));     /* Stationary Distribution */
      for(j=0; j<n_y; j++) {
	      ME->P_nml[i][j]  = (double *)calloc(n_z, sizeof(double));   
	      ME->PS_nml[i][j] = (double *)calloc(n_z, sizeof(double));   /* Stationary Distribution */
      }
    }
    
    ME->P_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->P_m_Marginal = (double *)calloc(n_y, sizeof(double));
    ME->P_l_Marginal = (double *)calloc(n_z, sizeof(double));

    ME->PS_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->PS_m_Marginal = (double *)calloc(n_y, sizeof(double));
    ME->PS_l_Marginal = (double *)calloc(n_z, sizeof(double));
  }
  else {
    printf(" This structure is only prepared to accept a maximum number of dimensions.\n");
    printf(" This number of defined in MODEL.h super header file\n");
    printf(" Your probability distribution has %d dimensions!!!\n", n_DIMENSION);
    printf(" Notice that when n_DIMENSION is bigger than 3, the probability distributions\n");
    printf(" are not defined as tensor objects, but as vectors over configurations,\n");
    printf(" where y[i] is the probability corresponding to the i-th configuration.\n"); 
    i = ME_n_DIMENSION_MAXIMUM;
    printf(" If n_DIMENSION (%d) is bigger than ME_n_DIMENSION_MAXIMUM (%d) the program will exit\n", 
    n_DIMENSION, i);
    assert(n_DIMENSION <= ME_n_DIMENSION_MAXIMUM);
  }

  ME->n_D = (int *)calloc(n_DIMENSION, sizeof(int));

  ME->Marginal_Probability_Label = (char **)calloc(n_DIMENSION, sizeof(char *) );
  for(i = 0; i<n_DIMENSION; i++) 
    ME->Marginal_Probability_Label[i] = (char *)calloc(100, sizeof(char));

  ME->Vector_Model_Variables = (double *)calloc(n_DIMENSION, sizeof(double)); 

  ME->TabConfi   = (int **)calloc(No_of_CONFIGURATIONAL_STATES, sizeof(int *));
  ME->Co = (configuration **)calloc(No_of_CONFIGURATIONAL_STATES, sizeof(configuration *));
  for(i = 0; i<No_of_CONFIGURATIONAL_STATES; i++) {
    ME->TabConfi[i] = (int *)calloc(n_DIMENSION, sizeof(int));
    ME->Co[i] = (configuration *)calloc(1, sizeof(configuration));
    ME->Co[i]->anDw = (int *)calloc(n_DIMENSION + 1, sizeof(int));
    ME->Co[i]->anUp = (int *)calloc(n_DIMENSION + 1, sizeof(int));  
    ME->Co[i]->n = (int *)calloc(n_DIMENSION, sizeof(int)); 
  }
}

void Master_Equation_Configurational_State_Setup (Master_Equation * ME )
{
  int i, D, n_0; 
  unsigned long long int S;

#if defined DIFFUSION_HII_nD
  configuration ** Co = ME->Co;
  int ** Configuration_Table = ME->TabConfi;

  Parameter_Table * Table = (Parameter_Table *)ME->Table;

  D   = ME->n_DIMENSION; 
  n_0 = Table->TOTAL_No_of_CONSUMERS; 

  S = General_Configuration_Table(D, n_0, Configuration_Table);

  assert( ME->No_of_CONFIGURATIONAL_STATES == (int)S );

  Generating_Network_of_Configurations(Co, Configuration_Table, D, S, n_0);
#endif
}

void Master_Equation_Free ( Master_Equation * ME )
{
  int i, j, k;

  Parameter_Table * Table = (Parameter_Table *)ME->Table; 
  
  int n_Time = Table->T->I_Time; 
  
  for(k=0; k<n_Time; k++) {  
    for(i=0; i < ME->n_DIMENSION; i++)
      free(ME->MPD_T[k][i]);
    
    free(ME->MPD_T[k]);
  }
  free(ME->MPD_T);  

  for(k=0; k<n_Time; k++) 
    free( ME->Evolving_Distribution[k] );
  free( ME->Evolving_Distribution );
 
  free( ME->Probability_Distribution );
  free( ME->Probability_Distribution_Time_0 );

  for(i=0; i<ME->n_DIMENSION; i++) {
      free(ME->MPD[i]);
      free(ME->MPD_S[i]);
  }
  free(ME->MPD);
  free(ME->MPD_S);
  
  if (ME->n_DIMENSION == 1) {
    free( ME->P_n );
    free(ME->P_n_Marginal);
    free(ME->PS_n_Marginal);
  }
  else if (ME->n_DIMENSION == 2) { 
    
    for(i=0; i<ME->n_x; i++) {
      free (ME->P_nm[i]);
      free (ME->PS_nm[i]);
    }
    free (ME->P_nm); free(ME->PS_nm); 

    free(ME->P_n_Marginal);
    free(ME->P_m_Marginal);

    free(ME->PS_n_Marginal);
    free(ME->PS_m_Marginal);
  }
  else if (ME->n_DIMENSION == 3) {
    
    for(i=0; i<ME->n_x; i++) {
      for(j=0; j<ME->n_y; j++) {
	      free (ME->P_nml[i][j]);
	      free (ME->PS_nml[i][j]);
      }
      free( ME->P_nml[i]);
      free( ME->PS_nml[i]); 
    }
    free (ME->P_nml); free(ME->PS_nml); 

    free(ME->P_n_Marginal);
    free(ME->P_m_Marginal);
    free(ME->P_l_Marginal);

    free(ME->PS_n_Marginal);
    free(ME->PS_m_Marginal);
    free(ME->PS_l_Marginal);
   
  }
  else {
    printf(" This structure is only prepared to accept a maximum number of dimensions.\n");
    printf(" This number is defined in MODEL.h super header file\n");
    printf(" Your probability distribution has %d dimensions!!!\n", ME->n_DIMENSION);
    printf(" Notice that when n_DIMENSION is bigger than 3, the probability distributions\n");
    printf(" are not defined as tensor objects, but as vectors over configurations,\n");
    printf(" where y[i] is the probability corresponding to the i-th configuration.\n"); 
    i = ME_n_DIMENSION_MAXIMUM;
    printf(" If n_DIMENSION (%d) is bigger than ME_n_DIMENSION_MAXIMUM (%d) the program will exit\n", 
    ME->n_DIMENSION, i);
    assert(ME->n_DIMENSION <= ME_n_DIMENSION_MAXIMUM);
  }

  free(ME->n_D);

  for(i = 0; i<ME->n_DIMENSION; i++) free(ME->Marginal_Probability_Label[i]);
  free(ME->Marginal_Probability_Label); 

  free(ME->Vector_Model_Variables); 

  for(i = 0; i<ME->No_of_CONFIGURATIONAL_STATES; i++) {
    free(ME->TabConfi[i]);
    free(ME->Co[i]->anDw);
    free(ME->Co[i]->anUp);
    free(ME->Co[i]->n);
    free(ME->Co[i]);
  }
  free(ME->TabConfi);
  free(ME->Co);

  free( ME );
}

void Master_Equation_Initialization (Master_Equation * ME,
				                            int No_of_CONFIGURATIONAL_STATES,
				                            int n_DIMENSION,
				                            int n_x,
				                            int n_y,
				                            int n_z) 
{
  int i;
  ME->No_of_CONFIGURATIONAL_STATES = No_of_CONFIGURATIONAL_STATES;
  ME->n_DIMENSION                  = n_DIMENSION;
  ME->n_x                          = n_x;
  ME->n_y                          = n_y;
  ME->n_z                          = n_z; 

  for(i=0; i<ME->n_DIMENSION; i++) { 
    if( i == 0 )  ME->n_D[i] = ME->n_x; 
    if( i == 1 )  ME->n_D[i] = ME->n_y; 
    if( i == 2 )  ME->n_D[i] = ME->n_z; 
    else          ME->n_D[i] = ME->n_x; 
    /* For instance, all dimension are equal to "1 + ME->Table->TOTAL_No_of_CONSUMERS"  
       for the case of the model: DIFFUSION_HII_nD
    */  
   }
}

/// @brief 
/// @param ME 
/// @param No_of_CONFIGURATIONAL_STATES 
/// @param n_DIMENSION 
/// @param n_x 
/// @param n_y 
/// @param n_z 
void Master_Equation_Stationary_Distribution_Allocation ( Master_Equation * ME,
							                                            int No_of_CONFIGURATIONAL_STATES,
							                                            int n_DIMENSION,
							                                            int n_x,
							                                            int n_y,
							                                            int n_z )
{
  /* This function should only be used where there is an interest in having a separete
     structure of Master_Equation type with the stationary probability distribution. 

     Input arguments:

     . n_DIMENSION, dimension of the mutivariate probability distribution. This is can 
     be maximum 3, this is, P(n, m, l) 
   
     . No_of_CONFIGURATIONAL_STATES, which is equal to (or less than ) the product of 
     n_x * n_y * n_z
     
     . n_x, n_y, n_z, are the maximum sizes for each dimension. 

     If the natural support of P(n, m, l) does not go from n=0 to n_x, m=0 to n_y, and 
     l=0 to n_z, then P(n, m, l) will take zero values out of its actual support (see 
     calloc() calls for zero allocations below).    
  */ 
  int i, j;
  Parameter_Table * Table = (Parameter_Table *)ME->Table; 
  
  assert(No_of_CONFIGURATIONAL_STATES < MAX_No_of_CONFIGURATIONAL_STATES);

  ME->Probability_Distribution=(double *)calloc(No_of_CONFIGURATIONAL_STATES, sizeof(double));

  ME->MPD_S = (double **)calloc(n_DIMENSION, sizeof(int *));
  for(i=0; i<n_DIMENSION; i++)
    ME->MPD_S[i] = (double *)calloc(Table->TOTAL_No_of_CONSUMERS+1, sizeof(int));

  if (n_DIMENSION == 1) {

    ME->PS_n = (double *)calloc(n_x, sizeof(double) );       /* Stationary Distribution */
    ME->PS_n_Marginal = (double *)calloc(n_x, sizeof(double));
  }
  else if (n_DIMENSION == 2) {

    ME->PS_nm = (double **)calloc(n_x, sizeof(double *));    /* Stationary Distribution */
    for(i=0; i<n_x; i++) {
      ME->PS_nm[i] = (double *)calloc(n_x, sizeof(double));  /* Stationary Distribution */
    }
  
    ME->PS_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->PS_m_Marginal = (double *)calloc(n_y, sizeof(double));
  }
  else if (n_DIMENSION == 3) {
    
    ME->PS_nml = (double ***)calloc(n_x, sizeof(double **));        /* Stationary Distribution */
    for(i=0; i<n_x; i++) {
      ME->PS_nml[i] = (double **)calloc(n_y, sizeof(double *));     /* Stationary Distribution */
      for(j=0; j<n_y; j++) {
	    ME->PS_nml[i][j] = (double *)calloc(n_z, sizeof(double));   /* Stationary Distribution */
      }
    }
    
    ME->PS_n_Marginal = (double *)calloc(n_x, sizeof(double));
    ME->PS_m_Marginal = (double *)calloc(n_y, sizeof(double));
    ME->PS_l_Marginal = (double *)calloc(n_z, sizeof(double));
  }
  else {
    printf(" This structure is only prepared to accept a maximum number of dimensions.\n");
    printf(" This number of defined in MODEL.h super header file\n");
    printf(" Your probability distribution has %d dimensions!!!\n", n_DIMENSION);
    printf(" Notice that when n_DIMENSION is bigger than 3, the probability distributions\n");
    printf(" are not defined as tensor objects!!!\n"); 
    assert(n_DIMENSION <= ME_n_DIMENSION_MAXIMUM);
  }

  ME->Marginal_Probability_Label = (char **)calloc(n_DIMENSION, sizeof(char *) );
  for(i = 0; i<n_DIMENSION; i++) 
    ME->Marginal_Probability_Label[i] = (char *)calloc(100, sizeof(char));

  ME->Vector_Model_Variables = (double *)calloc(n_DIMENSION, sizeof(double));  


}

void Master_Equation_Stationary_Distribution_Free ( Master_Equation * ME )
{
  int i, j;

  for(i=0; i<ME->n_DIMENSION; i++) 
    free(ME->MPD_S[i]);
  free(ME->MPD_S);

  free(ME->Probability_Distribution);

  if (ME->n_DIMENSION == 1) {
  
    free(ME->PS_n_Marginal);
  }
  else if (ME->n_DIMENSION == 2) { 
    
    for(i=0; i<ME->n_x; i++) {
  
      free (ME->PS_nm[i]);
    }
    free(ME->PS_nm); 

    free(ME->PS_n_Marginal);
    free(ME->PS_m_Marginal);
  }
  else if (ME->n_DIMENSION == 3) {
    
    for(i=0; i<ME->n_x; i++) {
      for(j=0; j<ME->n_y; j++) {
	
	free (ME->PS_nml[i][j]);
      }
      free( ME->PS_nml[i]); 
    }
    free(ME->PS_nml); 
  }
  else {
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", ME->n_DIMENSION); 
    assert(ME->n_DIMENSION < 4);
  }

  for(i = 0; i<ME->n_DIMENSION; i++) free(ME->Marginal_Probability_Label[i]);
  free(ME->Marginal_Probability_Label); 

  free(ME->Vector_Model_Variables); 
  
  free( ME );
}

void Stationary_Distribution_Allocation_Initialization( Parameter_Table * Table )
{
 
  int No_of_CONFIGURATIONAL_STATES;
  int n_DIMENSION;
  int n_x, n_y, n_z;

  Model_Parameters_Master_Equation(Table,
				   &No_of_CONFIGURATIONAL_STATES,
				   &n_DIMENSION,
				   &n_x, &n_y, &n_z);
  
  Master_Equation * MEq = (Master_Equation *)calloc( 1, sizeof(Master_Equation) );
  /* MEq and Table will be two data structures that point to each other */
  Table->MEq = MEq;           
  MEq->Table = Table;
  
  Master_Equation_Stationary_Distribution_Allocation ( MEq,
						       No_of_CONFIGURATIONAL_STATES,
						       n_DIMENSION,
						       n_x, n_y, n_z );
 
  Master_Equation_Initialization ( MEq,
				   No_of_CONFIGURATIONAL_STATES,
				   n_DIMENSION,
				   n_x, n_y, n_z );
  
  Master_Equation_Configurational_State_Setup ( MEq );

  Labels_for_Marginal_Probabilities( Table ); 
}

void Stationary_Distribution_Free( Parameter_Table * Table )
{
    Master_Equation_Stationary_Distribution_Free ( Table->MEq );
}
