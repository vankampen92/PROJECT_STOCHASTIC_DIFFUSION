#include <MODEL.h>

// #define EIGEN_VERBOSE  

void E_I_G_E_N___V_A_L_U_E___C_A_L_C_U_L_A_T_I_O_N ( double * y_Sol,  int K, int W, 
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
      evaluating_JACOBIAN_Matrix(m, y_Sol, 0., W, P);
      
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
