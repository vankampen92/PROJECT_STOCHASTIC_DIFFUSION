#include <MODEL.h>

#if defined VERBOSE
#define EIGEN_VERBOSE
#endif

void GSL_Eigenvalue_Calculation ( double * y_Sol,  int N, 
				  Parameter_Table * P, 
				  double * l_re, double * l_im )
{
  /* Input Arguments:
     
     - (l_re[i], l_im[i]) from i = 0 to i = N-1
     - y_Sol[i], from i=0 to i=N-1. System Resting Equilibrium 
     - The dimension of the system is N.
     - Parameter Table * P.

     Output parameters:
     
     - (l_re[i], l_im[i]) from i=1 to i=K+1. Set of complex-valued eigenvalues. 
     This function requires previous allocation of these arrays. 
  */
  int i; 
  int n;
  int K;

  K = N-1; /* Convention: Index of the latest stage */
  
  gsl_matrix * m = gsl_matrix_alloc(K+1, K+1); gsl_matrix_set_zero(m);
  n = K + 1;
  /* Calculating the Jacobian at the Fixed Point */
  JACOBIAN_Matrix(m, y_Sol, 0., K, P);

#if defined EIGEN_VERBOSE  
      printf(" Right after setting Jacobian Matrix()... "); // Print_Press_Key(1,0,"."); 
#endif         
      /*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* BEGIN : Controling that the input matrix has no GSL_NAN *  */    
      /*   END :                                                    */
 
#if defined EIGEN_VERBOSE  
      // printf("\nInitial Jacobian Matrix...\n");
      // show_a_view_FloatMatrix(mm, 1, K+1, 1, K+1); Print_Press_Key(1,0,".");
#endif
      
      gsl_matrix_view M = gsl_matrix_view_array (m->data, n, n);

      gsl_vector_complex * eval = gsl_vector_complex_alloc (n);

      gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (n);
      
      /* Balancing is performed */
      gsl_eigen_nonsymm_params (0, 1, w);

      gsl_eigen_nonsymm (&M.matrix, eval, w);

      gsl_eigen_nonsymm_free (w);

#if defined EIGEN_VERBOSE  
      printf(" Eigenvalues are: \n");
      for (i = 0; i < n; i++) {

             gsl_complex eval_i = gsl_vector_complex_get (eval, i);
             
             printf ("eigenvalue = %g + i %g\n",
                     GSL_REAL(eval_i), GSL_IMAG(eval_i));
             
      }
      // Print_Press_Key(1,0,".");     
#endif 
      
      /* BEGIN : Generating output float vectors of full dimension */
      for( i = 0; i < N; i++ ){

       	  gsl_complex eval_i = gsl_vector_complex_get (eval, i);

	  l_re[i] = GSL_REAL(eval_i);
	  l_im[i] = GSL_IMAG(eval_i);
      } 
      /*  END : --------------------------------------------------*/


#if defined EIGEN_VERBOSE  	  
      // printf("Eigen Vectors ...\n");
      // int showing = showing_eigenValues(l_re, l_im, K+1); Print_Press_Key(1,0,".");
#endif
      
      /* fp = fopen("eigenValues.dat", "w"); */
      /* 	int saving = saving_eigenValues(fp, l_re, l_im, K+1); */
      /* 	fclose(fp); */
    
      
      gsl_vector_complex_free(eval);       
      gsl_matrix_free(m);
}

void Dominant_Eigenvalue_Calculation(double * Y1, double * Y2, int N,
				     int * Index_Value_D, int * Index_Value_S)
{
  /* This function takes complex-valued eigenvalues of the form:
     
                                   Y1[k] + i Y2[k] 
                                   from k=0 to N-1
     
     and search for the dominant eigenvalue and the subdominant eigen value. The out
     put is given in terms of pointers to the indeces that allow to retrieve the dominant 
     and the subdominant eigenvalues: 

     . Index_Value_D
     . Index_Value_S
  */ 

  int i, k, n;
  double Real_V, Imm_V;

  double * S_Y1 = (double *)calloc(N-1, sizeof(double) );
  double * S_Y2 = (double *)calloc(N-1, sizeof(double) ); 
  int * Index_S = (int *)calloc(N-1, sizeof(int));
  
  Real_V = Y1[0]; Imm_V = Y2[0]; k = 0; 
  for (i = 0; i < N; i++) {
    
    printf ("EigenValue[%d] = %g + i %g\n", i, Y1[i], Y2[i]);
    
    if( Real_V < Y1[i] ) {
      Real_V = Y1[i];
      Imm_V  = Y2[i]; 
      k = i;
    }
  }
  * Index_Value_D = k;
  printf("Valor Propi Dominant: V[%d] = = %g + i %g\n", k, Y1[k], Y2[k]);
  
  n=0; 
  for (i = 0; i < N; i++) {
    if (i != k) {
      S_Y1[n] = Y1[i];
      S_Y2[n] = Y2[i];
      Index_S[n] =  i; 
      n++; 
    }
  }
  
  Real_V = S_Y1[0]; Imm_V = S_Y2[0]; k = 0;  
  for (i = 0; i < n; i++) {
  
    if( Real_V < S_Y1[i] ) {
	  Real_V = S_Y1[i];
	  Imm_V  = S_Y2[i]; 
	  k = i;
    }
  }
  * Index_Value_S = Index_S[k];
      
  printf("Valor Propi Sub-Dominant: V[%d] = %g + i %g\n", Index_S[k], S_Y1[k], S_Y2[k]);

  k = * Index_Value_S; 
  printf("Valor Propi Sub-Dominant: V[%d] = %g + i %g\n", k, Y1[k], Y2[k]);
  
  //Print_Press_Key(1,0,".");  
  
  free(S_Y1); free(S_Y2); free(Index_S); 
}  


/* void E_I_G_E_N___V_A_L_U_E___C_A_L_C_U_L_A_T_I_O_N ( double * y_Sol,  int K, int W,  */
/* 						     Parameter_Table * P,  */
/* 						     float * l_re, float * l_im ) */
/* { */
/*   /\* Input Arguments: */
     
/*      - (l_re[i], l_im[i]) from i = 1 to i = K+1 */
/*      - y_Sol[i], from i=0 to i=K. System Resting Equilibrium  */
/*      - The dimension of the system is K + 1 (W + 1). */
/*      - Parameter Table * P. */

/*      Output parameters: */
     
/*      - (l_re[i], l_im[i]) from i=1 to i=K+1. Set of complex-valued eigenvalues.  */
/*      This function requires previous allocation of these arrays.  */
/*   *\/ */
/*   int i;  */
/*   int n; */

/*   gsl_matrix * m = gsl_matrix_alloc(K+1, K+1); gsl_matrix_set_zero(m); */
/*   n = K + 1; */
/*   /\* Calculating the Jacobian at the Fixed Point *\/ */
/*   JACOBIAN_Matrix(m, y_Sol, 0., W, P); */

/* #if defined EIGEN_VERBOSE   */
/*       printf(" Right after setting Jacobian Matrix()... "); Print_Press_Key(1,0,".");  */
/* #endif          */
/*       /\*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\/ */
/*       /\* BEGIN : Controling that the input matrix has no GSL_NAN *  *\/     */
/*       /\*   END :                                                    *\/ */
 
/* #if defined EIGEN_VERBOSE   */
/*       // printf("\nInitial Jacobian Matrix...\n"); */
/*       // show_a_view_FloatMatrix(mm, 1, K+1, 1, K+1); Print_Press_Key(1,0,"."); */
/* #endif */
      
/*       gsl_matrix_view M = gsl_matrix_view_array (m->data, n, n); */

/*       gsl_vector_complex * eval = gsl_vector_complex_alloc (n); */

/*       gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (n); */
      
/*       /\* Balancing is performed *\/ */
/*       gsl_eigen_nonsymm_params (0, 1, w); */

/*       gsl_eigen_nonsymm (&M.matrix, eval, w); */

/*       gsl_eigen_nonsymm_free (w); */

/* #if defined EIGEN_VERBOSE   */
/*       printf(" Eigenvalues are: \n"); */
/*       for (i = 0; i < n; i++) { */

/*              gsl_complex eval_i = gsl_vector_complex_get (eval, i); */
             
/*              printf ("eigenvalue = %g + i %g\n", */
/*                      GSL_REAL(eval_i), GSL_IMAG(eval_i)); */
             
/*       } */
/*       Print_Press_Key(1,0,".");      */
/* #endif  */
      
/*       /\* BEGIN : Generating output float vectors of full dimension *\/ */
/*       set_to_value_float_1_N(l_re, K+1, 0.0);  */
/*       set_to_value_float_1_N(l_im, K+1, 0.0);  */

/*       for( i = 1; i<= (K+1); i++ ){ */

/*        	  gsl_complex eval_i = gsl_vector_complex_get (eval, i-1); */

/* 	  l_re[i] = (float)GSL_REAL(eval_i); */
/* 	  l_im[i] = (float)GSL_IMAG(eval_i); */
/*       }  */
/*       /\*  END : --------------------------------------------------*\/ */


/* #if defined EIGEN_VERBOSE  	   */
/*       // printf("Eigen Vectors ...\n"); */
/*       // int showing = showing_eigenValues(l_re, l_im, K+1); Print_Press_Key(1,0,"."); */
/* #endif */
      
/*       /\* fp = fopen("eigenValues.dat", "w"); *\/ */
/*       /\* 	int saving = saving_eigenValues(fp, l_re, l_im, K+1); *\/ */
/*       /\* 	fclose(fp); *\/ */
    
      
/*       gsl_vector_complex_free(eval);        */
/*       gsl_matrix_free(m); */
/* } */

