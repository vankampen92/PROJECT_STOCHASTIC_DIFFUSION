double x, K_R;

  K_R = (double)Table->K_R;

  /* F_R( y; parameters) */

  x= -Table->Lambda_R_0 -Table->Delta_R_0 +Table->Beta_R -2.0*Table->Beta_R*y[R]/K_R -Table->Alpha_C_0 *y[A]/K_R;
  gsl_matrix_set(m, R, 0, x);

  x= -Table->Alpha_C_0 *y[R]/K_R;
  gsl_matrix_set(m, R, 1, x);

  x= 0.0;
  gsl_matrix_set(m, R, 2, x);

  x= 0.0;
  gsl_matrix_set(m, R, 3, x);

  /* F_A( y; parameters) */

  x= -Table->Alpha_C_0 *y[A]/K_R;
  gsl_matrix_set(m, A, 0, x);

  x= -Table->Delta_C_0 -Table->Alpha_C_0 *y[R]/K_R -Table->Chi_C_0 *y[RA]/K_R;
  gsl_matrix_set(m, A, 1, x);

  x= +2.0*Table->Nu_C_0 -Table->Chi_C_0 *y[A]/K_R;
  gsl_matrix_set(m, A, 2, x);

  x= +Table->Eta_C_0;
  gsl_matrix_set(m, A, 3, x);

  /* F_AR( y; parameters) */
  x= Table->Alpha_C_0 *y[A]/K_R;
  gsl_matrix_set(m, RA, 0, x);

  x= Table->Alpha_C_0 *y[R]/K_R -Table->Chi_C_0 *y[RA]/K_R;  
  gsl_matrix_set(m, RA, 1, x);

  x= -Table->Nu_C_0 -Table->Chi_C_0 *y[A]/K_R;
  gsl_matrix_set(m, RA, 2, x);

  x= +Table->Eta_C_0;
  gsl_matrix_set(m, RA, 3, x);

/* F_ARA( y; parameters) */

  x= 0.0;
  gsl_matrix_set(m, ARA, 0, x);

  x= Table->Chi_C_0 *y[RA]/K_R;
  gsl_matrix_set(m, ARA, 1, x);

  x= Table->Chi_C_0 *y[A]/K_R;
  gsl_matrix_set(m, ARA, 2, x);

  x= -Table->Eta_C_0;
  gsl_matrix_set(m, ARA, 3, x);

  if( Table->No_of_CELLS > 1 ) {
    /* This Jacobian Matrix evaluation only works 
       for single-patch systems 
    */
    assert( Table->No_of_CELLS == 1 ); 
  }  

