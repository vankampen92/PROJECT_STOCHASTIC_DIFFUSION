#include <MODEL.h>

#define EIGEN_VERBOSE  
// #define DYNAMIC___VERBOSE

void NR_Eigenvalue_Calculation ( double * y_Sol,  int N, 
				 Parameter_Table * P, 
				 double * l_re, double * l_im )
{
  int i, K;
  K  = N-1;

  float * l_re_f = (float *)calloc( N, sizeof(float) );
  float * l_im_f = (float *)calloc( N, sizeof(float) ); 

  NR_Eigenvalue_Calculation_float ( y_Sol,  K, 0, 
				    P, 
				    l_re_f, l_im_f );
  for(i = 0; i<N; i++) {
    l_re[i] = (double)l_re_f[i+1];
    l_im[i] = (double)l_im_f[i+1];
  }

  free(l_re_f); free(l_im_f); 
}

void NR_Eigenvalue_Calculation_float ( double * y_Sol,  int K, int W, 
				       Parameter_Table * P, 
				       float * l_re, float * l_im )
{
  gsl_matrix * m;
  float **mm; 
  
  float nr___x_min = P->nr___x_min;
  float nr___x_MAX = P->nr___x_MAX;
     
  m     = gsl_matrix_alloc(K+1, K+1); gsl_matrix_set_zero(m);
  mm = matrix(1,K+1, 1,K+1);          set_matrix_to_value_float_1_Nx_1_Ny (mm, K+1, K+1, 0.0); 
	     
  set_to_value_float_1_N(l_re, K+1, 0.0); 
  set_to_value_float_1_N(l_im, K+1, 0.0); 
	
     /* Calculating the Jacobian at the Fixed Point */
      JACOBIAN_Matrix(m, y_Sol, 0., W, P);
      
      int gsl_to_NR = gsl_matrix_to_NR_matrix(m, mm, K+1, K+1);
      
      /*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* BEGIN : Controling that the input matrix has no GSL_NAN *  */    
      /*   END :                                                    */
      
#if defined EIGEN_VERBOSE  
      printf("\nInitial Jacobian Matrix...\n");
      if ( Well_Defined_Matrix_Elements(mm, K+1, nr___x_min, nr___x_MAX) == 1.0 ){
	printf(" All matrix elements are well-defined after\n");
      }
      else {
	printf(" Further calculations will have NaN problems!!!\n");
	printf(" Some matrix elements does not seem to be well-defined after\n");
      }
      
      printf(" the Jacobian evaluation by using :\n");
      printf(" evaluating_JACOBIAN_Matrix(m, y_Sol->data, 0., W, P)\n");
      show_a_view_FloatMatrix(mm, 1, K+1, 1, K+1); Print_Press_Key(1,0,".");
      Print_Press_Key(1,0,".");
#endif
      
      balanc(mm, K+1); 
      
#if defined EIGEN_VERBOSE
      printf("Balaced Jacobian Matrix...\n");
      if ( Well_Defined_Matrix_Elements(mm, K+1, nr___x_min, nr___x_MAX) == 1.0 ){
	printf(" All matrix elements are well-defined after\n");
      }
      else {
	printf(" Further calculations will have NaN problems!!!\n");
	printf(" Some matrix elements does not seem to be well-defined after\n");
      }
      printf(" balancing the Jacobian by using :\n");
      printf(" balanc(mm, K+1) \n");
      show_a_view_FloatMatrix(mm, 1, K+1, 1, K+1); Print_Press_Key(1,0,".");
      Print_Press_Key(1,0,".");
#endif
      
      elmhes(mm, K+1);
      
#if defined EIGEN_VERBOSE
      printf("Corresponding Hessenberg Form of the Jacobian Matrix...\n");
      if ( Well_Defined_Matrix_Elements(mm, K+1, nr___x_min, nr___x_MAX) == 1.0 ){
	printf(" All matrix elements are well-defined after\n");
      }
      else {
	printf(" Further calculations will have NaN problems!!!\n");
	printf(" Some matrix elements does not seem to be well-defined after\n");
      }
      printf(" calculating the Hessenberg matrix from the balanced Jacobian by using :\n");
      printf(" elmhes(mm, K+1) \n");
      show_a_view_FloatMatrix(mm, 1, K+1, 1, K+1); Print_Press_Key(1,0,".");
      Print_Press_Key(1,0,".");
#endif

      hqr(mm, K+1, l_re, l_im);
      
#if defined EIGEN_VERBOSE  	  
      printf("Eigen Vectors ...\n");
      int showing = showing_eigenValues(l_re, l_im, K+1); Print_Press_Key(1,0,".");
#endif
      
      /* fp = fopen("eigenValues.dat", "w"); */
      /* 	int saving = saving_eigenValues(fp, l_re, l_im, K+1); */
      /* 	fclose(fp); */
    
      free_matrix(mm, 1,K+1, 1,K+1);
      gsl_matrix_free(m);
}

int showing_eigenValues(float *l_re, float *l_im, int n)
{
  int i;
  
  for(i=1; i<=n; i++){
    printf("( %f ) + i( %f )\n", l_re[i], l_im[i]);
  }
  
  return(0);
}

int fprintf_to_File_Matrix_gsl(FILE * Out, gsl_matrix * A, int MM, int NN)
{
  int i,j;
  for (i=0;i<MM;i++){
    for (j=0;j<NN;j++){
      fprintf(Out,"%g\t",gsl_matrix_get(A,i,j));
    }
    fprintf(Out,"\n");
  }
  fprintf(Out,"\n");
  return 0;
}

int gsl_matrix_to_NR_matrix(gsl_matrix * A, float **a,  int MM, int NN)
{
  /* Utility to transform a gsl matrix into a numerical recipes float matrix */
  int i,j;
  for (i=1;i<=MM;i++){
    for (j=1;j<=NN;j++){
      a[i][j] = (float)gsl_matrix_get(A,i-1,j-1);
    }
  }
  return 0;
}

double Well_Defined_Matrix_Elements(float **mm, int N, float xmin, float xmax)
{
  /* 
     xmin and xmax are the largest and tiniest float representation
     of a number the machine can handle without producing overflow.
     Out of this range overflow is returned. 
     
     This function assesses whether all elements in matrix mm
     fall within this range and therefore they are all number
     representations the machine will be able to handle without 
     giving neither further problems nor surprises.
  */
  int i,j, count;
  double x;
  
  x = 0.0; 
  count = 0;
  
  for(i=1; i<=N; i++){
    for(j=1; j<=N; j++){
      
      if ( mm[i][j] != 0.0 ) {
	//if(mm[i][j] < xmax && fabsf( mm[i][j] ) >= xmin ){
	if(mm[i][j] < xmax && fabsf( mm[i][j] ) >= 0.0 ){
	  count++;	
#if defined DYNAMIC___VERBOSE
	  printf("Matrix element: A(%d, %d) = %g\n", i,j, mm[i][j]);
#endif
	}
	else{
#if defined DYNAMIC___VERBOSE
	  printf("Ill-defiend matrix element: A(%d, %d) = %g\n", i,j, mm[i][j]);
#endif
	}
      }
      else {
	count++;
#if defined DYNAMIC___VERBOSE
	printf("Nul matrix element: A(%d, %d) = %g\n", i,j, mm[i][j]);
#endif
      }
    }						
  }
  
  if(count == N * N) 
    x = 1.0;
  else
    x = 0.0;

#if defined DYNAMIC___VERBOSE 
  printf("count = %d\t Matrix dimension (N) = %d\t N * N = %d\n", count, N, N*N);
#endif

  return(x);
}
