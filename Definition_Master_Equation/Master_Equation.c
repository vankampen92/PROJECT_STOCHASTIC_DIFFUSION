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
     calloc() calls for zero allocations below).    
  */ 
  int i, j;
  Parameter_Table * Table = (Parameter_Table *)ME->Table; 
  
  assert(No_of_CONFIGURATIONAL_STATES < MAX_No_of_CONFIGURATIONAL_STATES);
  
  ME->Probability_Distribution = (double *)calloc(No_of_CONFIGURATIONAL_STATES, sizeof(double) );

  ME->Probability_Distribution_Time_0 = (double *)calloc(No_of_CONFIGURATIONAL_STATES,
							sizeof(double));

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
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", n_DIMENSION); 
    assert(n_DIMENSION < 4);
  }

  ME->Marginal_Probability_Label = (char **)calloc(n_DIMENSION, sizeof(char *) );
  for(i = 0; i<n_DIMENSION; i++) 
    ME->Marginal_Probability_Label[i] = (char *)calloc(100, sizeof(char));

  ME->Vector_Model_Variables = (double *)calloc(n_DIMENSION, sizeof(double)); 
  
}

void Master_Equation_Free ( Master_Equation * ME )
{
  int i, j;

  free( ME->Probability_Distribution );
  free( ME->Probability_Distribution_Time_0 );
  
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
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", ME->n_DIMENSION); 
    assert(ME->n_DIMENSION < 4);
  }

  for(i = 0; i<ME->n_DIMENSION; i++) free(ME->Marginal_Probability_Label[i]);
  free(ME->Marginal_Probability_Label); 

  free(ME->Vector_Model_Variables); 
  
  free( ME );
}

void Master_Equation_Initialization (Master_Equation * ME,
				     int No_of_CONFIGURATIONAL_STATES,
				     int n_DIMENSION,
				     int n_x,
				     int n_y,
				     int n_z) 
{
  ME->No_of_CONFIGURATIONAL_STATES = No_of_CONFIGURATIONAL_STATES;
  ME->n_DIMENSION                  = n_DIMENSION;
  ME->n_x                          = n_x;
  ME->n_y                          = n_y;
  ME->n_z                          = n_z; 
}

void Master_Equation_Stationary_Distribution_Allocation ( Master_Equation * ME,
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
     calloc() calls for zero allocations below).    
  */ 
  int i, j;
  Parameter_Table * Table = (Parameter_Table *)ME->Table; 
  
  assert(No_of_CONFIGURATIONAL_STATES < MAX_No_of_CONFIGURATIONAL_STATES);
  
  
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
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", n_DIMENSION); 
    assert(n_DIMENSION < 4);
  }

  ME->Marginal_Probability_Label = (char **)calloc(n_DIMENSION, sizeof(char *) );
  for(i = 0; i<n_DIMENSION; i++) 
    ME->Marginal_Probability_Label[i] = (char *)calloc(100, sizeof(char));

  ME->Vector_Model_Variables = (double *)calloc(n_DIMENSION, sizeof(double));  
}

void Master_Equation_Stationary_Distribution_Free ( Master_Equation * ME )
{
  int i, j;

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
  
  Labels_for_Marginal_Probabilities( Table ); 
}

void Stationary_Distribution_Free( Parameter_Table * Table )
{
    Master_Equation_Stationary_Distribution_Free ( Table->MEq );
}
