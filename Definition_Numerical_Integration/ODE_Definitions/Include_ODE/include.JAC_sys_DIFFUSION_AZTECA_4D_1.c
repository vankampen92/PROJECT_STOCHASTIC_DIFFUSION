double x;

/* K_Q: Total Max No of Potential Nesting Trees (per local patch) */ 
double K_Q = (double)Table->K_R; 

  /* F_W( y; parameters) */
  x = -Table->Delta_R_0  -Table->Alpha_C_0* y[F]/K_Q; 
  gsl_matrix_set(m, W, 0, x);
 
  x = Table->Beta_R;
  gsl_matrix_set(m, W, 1, x);

  x= Table->Alpha_C_0 * y[W]/K_Q;
  gsl_matrix_set(m, W, 2, x);

  x= 0.0;
  gsl_matrix_set(m, W, 3, x);

  /* F_Q( y; parameters) */
  x = Table->Eta_R *(K_Q - y[Q])/K_Q; 
  gsl_matrix_set(m, Q, 0, x);

  x= -Table->Eta_R *y[W]/K_Q -Table->Delta_R_1;
  gsl_matrix_set(m, Q, 1, x);

  x= 0.0;
  gsl_matrix_set(m, Q, 2, x);

  x= 0.0;
  gsl_matrix_set(m, Q, 3, x);

  /* F_F( y; parameters) */
  x= 0.0;
  gsl_matrix_set(m, F, 0, x);

  x= 0.0;
  gsl_matrix_set(m, F, 1, x);

  x= -Table->Delta_C_0;
  gsl_matrix_set(m, F, 2, x);

  x= +Table->Nu_C_0;
  gsl_matrix_set(m, F, 3, x);

  /* F_WF( y; parameters) */
  x= Table->Alpha_C_0 *y[F]/K_Q;
  gsl_matrix_set(m, WF, 0, x);

  x= 0.0;
  gsl_matrix_set(m, WF, 1, x);

  x= Table->Alpha_C_0 *y[W]/K_Q ;  
  gsl_matrix_set(m, WF, 2, x);

  x= -Table->Nu_C_0 -Table->Delta_C_1;
  gsl_matrix_set(m, WF, 3, x);

  if( Table->No_of_CELLS > 1 ) {
    /* This Jacobian Matrix evaluation only works 
       for single-patch systems 
    */
    assert( Table->No_of_CELLS == 1 ); 
  }  

