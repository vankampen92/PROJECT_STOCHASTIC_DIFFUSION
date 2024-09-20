double x, K_R;

  K_R = (double)Table->K_R;

  /* F_RP( y; parameters) */
  x = -Table->Delta_R_1 - Table->Eta_R* (1.0 - y[R]/K_R); 
  gsl_matrix_set(m, RP, 0, x);
 
  x = Table->Beta_R +  Table->Eta_R * y[RP]/K_R;
  gsl_matrix_set(m, RP, 1, x);


  /* F_R( y; parameters) */
  x = Table->Eta_R *(1.0 - y[R]/K_R); 
  gsl_matrix_set(m, R, 0, x);

  x= -Table->Delta_R_0 -Table->Eta_R *y[RP]/K_R -Table->Alpha_C_0 *y[A]/K_R;
  gsl_matrix_set(m, R, 1, x);


  if( Table->No_of_CELLS > 1 || Table->No_of_RESOURCES > 1) {
    /* This Jacobian Matrix evaluation only works 
       for single-patch systems 
    */
    assert( Table->No_of_CELLS == 1 ); 
    assert( Table->No_of_RESOURCES == 1);
  }  

