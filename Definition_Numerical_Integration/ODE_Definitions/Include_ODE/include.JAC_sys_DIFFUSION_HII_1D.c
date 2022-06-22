double x, K_R, y_R;

  K_R = (double)Table->K_R;
  y_R = Table->TOTAL_No_of_RESOURCES; /* Resource density (in number of resource units) */

  /* F_A( y; parameters) */

  x = -Table->Alpha_C_0 *y_R/K_R - Table->Nu_C_0;
  gsl_matrix_set(m, A, 0, x);

  if( Table->No_of_CELLS > 1 ) {
    /* This Jacobian Matrix evaluation only works 
       for single-patch systems. The program 
       will stop here!!! 
    */
    printf("This Jacobian Matrix evaluation only works\n");
    printf("for single-patch systems. The program will\n");
    printf("stop here!!!\n");
    
    assert( Table->No_of_CELLS == 1 ); 
  }  

 




 




