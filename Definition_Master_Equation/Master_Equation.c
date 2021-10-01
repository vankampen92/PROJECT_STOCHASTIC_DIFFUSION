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
				  int n_z)
{
  /* Input arguments: 
     
     . n_DIMENSION, dimension of the mutivariate probability distribution. This is can 
     be maximum 3, this is, P(n, m, l) 
   
     . No_of_CONFIGURATIONAL_STATES, which is equal to the product of n_x * n_y * n_z
     
     . n_x, n_y, n_z, size for each dimension. 
     
  */ 
  int i, j;

  assert(No_of_CONFIGURATIONAL_STATES < MAX_No_of_CONFIGURATIONAL_STATES);

  assert( n_x * n_y * n_z == No_of_CONFIGURATIONAL_STATES );  
  
  ME->Probability_Distribution = (double *)caloc(No_of_CONFIGURATIONAL_STATES, sizeof(double) );

  ME->Probability_Distribution_Time_0 = (double *)caloc(No_of_CONFIGURATIONAL_STATES,
							sizeof(double));

  if (n_DIMENSION == 1)
    ME->P_n = (double *)calloc(n_x, sizeof(double) );

  else if (n_DIMENSION == 2) { 
    ME->P_nm = (double **)calloc(n_x, sizeof(double *)); 
    for(i=0; i<n_x; i++) 
      ME->P_nm[i] = (double *)calloc(n_y, sizeof(double)); 
  }
  else if (n_DIMENSION == 3) {
    ME->P_nml = (double ***)calloc(n_x, sizeof(double **)); 
    for(i=0; i<n_x; i++) {
      ME->P_nml[i] = (double **)calloc(n_y, sizeof(double *));
      for(j=0; j<n_y; j++) 
	ME->P_nml[i][j] = (double *)calloc(n_z, sizeof(double));
    }
  }
  else {
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", n_DIMENSION); 
    assert(n_DIMENSION < 4);
  }
  
}

void Master_Equation_Free ( Master_Equation * ME )
{
  int i, j;

  free( ME->Probability_Distribution );
  free( ME->Probability_Distribution_Time_0 );
  
  if (ME->n_DIMENSION == 1)
    free( ME->P_n );

  else if (ME->n_DIMENSION == 2) { 
    
    for(i=0; i<ME->n_x; i++) free (ME->P_nm[i]); 
    free (ME->P_nm);
  }
  else if (ME->n_DIMENSION == 3) {
    
    for(i=0; i<ME->n_x; i++) {
      for(j=0; j<ME->n_y; j++) free (ME->P_nml[i][j]); 
      free( ME->P_nml[i]); 
    }
    free (ME->P_nml);
  }
  else {
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", ME->n_DIMENSION); 
    assert(ME->n_DIMENSION < 4);
  }
  
  free( ME );
}

void Master_Equation_Initialization (Master_Equation * ME,
				     int No_of_CONFIGURATIONAL_STATES,
				     int n_DIMENSION,
				     int n_x,
				     int n_y,
				     int n_z) )
{
  ME->No_of_CONFIGURATIONAL_STATES = No_of_CONFIGURATIONAL_STATES;
  ME->n_DIMENSION                  = n_DIMENSION;
  ME->n_x                          = n_x;
  ME->n_y                          = n_y;
  ME->n_z                          = n_z; 
}

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME )
{
  int i;
  int n,m,l;
  int s; 

  double * y = ME->Probability_Distribution; 

  if (ME->n_DIMENSION == 1) {
    for(i=0; i<ME->n_x; i++) ME->P_n[i] = y[i];
  }  
  else if (ME->n_DIMENSION == 2) {
    
    for(i=0; i<ME->No_of_CONFIGURATION_STATES; i++) {
      n = i/ME->n_y;
      m = i%ME->n_y;
      ME->P_nm[n][m] = y[i];
    }
    
  }
  else if (ME->n_DIMENSION == 3) {

    for(i=0; i<ME->No_of_CONFIGURATION_STATES; i++) {
      n = i/(ME->n_y*ME->n_z);
      s = i%(ME->n_y*ME->n_z);
      m = s/ME->n_z;
      l = s%ME->n_z
	
      ME->P_nml[n][m][l] = y[i];
    }
    
  }
  else {
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", ME->n_DIMENSION); 
    assert(ME->n_DIMENSION < 4);
  }
}
