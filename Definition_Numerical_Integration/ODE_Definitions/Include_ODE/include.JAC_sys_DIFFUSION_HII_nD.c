double x;

  /* F_ik( y; parameters) */

  for(i=0; i<Table->No_of_RESOURCES; i++) {
    for(k=0; k<Table->No_of_RESOURCES; k++) {
      if( i == k ) {
        x = -(Table->Nu_Consumers[i] + Table->Theta_Consumers[k]);  
        gsl_matrix_set(m, i, k, x);
      }
      else {
        x = -Table->Theta_Consumers[k];  
        gsl_matrix_set(m, i, k, x);
      }
    }
  }

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

 




 




