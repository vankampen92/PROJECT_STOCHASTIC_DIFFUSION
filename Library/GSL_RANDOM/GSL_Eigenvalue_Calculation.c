#include "GSL_stat.h"

#if defined VERBOSE
#define EIGEN_VERBOSE
#endif 

void E_I_G_E_N___V_A_L_U_E___C_A_L_C_U_L_A_T_I_O_N ( double * y_Sol,  int K, int W, 
						     Parameter_Table * P, 
						     float * l_re, float * l_im )
{
  /* Input Arguments:
     
     - (l_re[i], l_im[i]) from i = 1 to i = K+1
     - y_Sol[i], from i=0 to i=K. System Resting Equilibrium 
     - The dimension of the system is K + 1 (W + 1).
     - Parameter Table * P.

     Output parameters:
     
     - (l_re[i], l_im[i]) from i=1 to i=K+1. Set of complex-valued eigenvalues. 
     This function requires previous allocation of these arrays. 
  */
  int i; 
  int n;

  gsl_matrix * m = gsl_matrix_alloc(K+1, K+1); gsl_matrix_set_zero(m);
  n = K + 1;
  /* Calculating the Jacobian at the Fixed Point */
  evaluating_JACOBIAN_Matrix(m, y_Sol, 0., W, P);

#if defined EIGEN_VERBOSE  
      printf(" Right after setting Jacobian Matrix()... "); Press_Key(); 
#endif         
      /*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* BEGIN : Controling that the input matrix has no GSL_NAN *  */    
      /*   END :                                                    */
 
#if defined EIGEN_VERBOSE  
      // printf("\nInitial Jacobian Matrix...\n");
      // show_a_view_FloatMatrix(mm, 1, K+1, 1, K+1); Press_Key();
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
      Press_Key();     
#endif 
      
      /* BEGIN : Generating output float vectors of full dimension */
      set_to_value_float_1_N(l_re, K+1, 0.0); 
      set_to_value_float_1_N(l_im, K+1, 0.0); 

      for( i = 1; i<= (K+1); i++ ){

       	  gsl_complex eval_i = gsl_vector_complex_get (eval, i-1);

	  l_re[i] = (float)GSL_REAL(eval_i);
	  l_im[i] = (float)GSL_IMAG(eval_i);
      } 
      /*  END : --------------------------------------------------*/


#if defined EIGEN_VERBOSE  	  
      // printf("Eigen Vectors ...\n");
      // int showing = showing_eigenValues(l_re, l_im, K+1); Press_Key();
#endif
      
      /* fp = fopen("eigenValues.dat", "w"); */
      /* 	int saving = saving_eigenValues(fp, l_re, l_im, K+1); */
      /* 	fclose(fp); */
    
      
      gsl_vector_complex_free(eval);       
      gsl_matrix_free(m);
}
