#include <MODEL.h>

double R0_Function_F( Parameter_Table * Table )
{
  /* Invasion criteria for flies into a population of ants, 
     when grows stably as a population 
  */
  int i;
  double b_W, e_Q, n_WF, n_F, a_W, a_F, d_W, R_W, R_F;
  double x, xW, xQ, xF, xWF; 
                                               /* Delta_W is  Delta_R_0 */
  b_W = Table->Beta_R/Table->Delta_R_0;        /* Delta_Q is  Delta_R_1 */
  e_Q = Table->Eta_R/Table->Delta_R_1;         /* Delta_F is  Delta_C_0 */
  n_WF = Table->Nu_C_0/Table->Delta_C_1;       /* Delta_WF is  Delta_C_1 */
  n_F  = Table->Nu_C_0/Table->Delta_C_0; 
  a_W  = Table->Alpha_C_0/Table->Delta_R_0;
  a_F  = Table->Alpha_C_0/Table->Delta_C_0;
  d_W  = Table->Delta_C_1/Table->Delta_R_0;

  assert_right_model_definition( Table );

  R_W = b_W * e_Q;                           /* Threshold for Ants Persistence */

  #include <Model_Variables_Code.Include.c>

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {
 
    /* Under these conditions, four Dynamic Regimes at Stationarity are possible */
    if ( R_W < 1.0 ){
      
      R_F = 0;
     
      printf("Self-maintainance of the ant population is not possible\n");
      printf("Therefore, no flies can enter!!!\t R_0 = %g", R_F);
    }
    else { 
      
      R_F = R_W / (1.0 + (e_Q /a_F)*(1.0 + 1.0/n_WF));
  }
  else {
    printf("R_0 is only possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, here, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g)\n", 
          Table->Lambda_C_0, Table->Lambda_R_0);
    printf("or the total number of patches is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }

  return(R_F);
}

double R0_Function_W( Parameter_Table * Table )
{
  /* Invasion criteria for ants into an empty landscape 
  */
  int i;
  double b_W, e_Q, n_WF, n_F, a_W, a_F, d_W, R_W, R_F;
  double x, xW, xQ, xF, xWF; 
                                               /* Delta_W is  Delta_R_0 */
  b_W = Table->Beta_R/Table->Delta_R_0;        /* Delta_Q is  Delta_R_1 */
  e_Q = Table->Eta_R/Table->Delta_R_1;         /* Delta_F is  Delta_C_0 */
  n_WF = Table->Nu_C_0/Table->Delta_C_1;       /* Delta_WF is  Delta_C_1 */
  n_F  = Table->Nu_C_0/Table->Delta_C_0; 
  a_W  = Table->Alpha_C_0/Table->Delta_R_0;
  a_F  = Table->Alpha_C_0/Table->Delta_C_0;
  d_W  = Table->Delta_C_1/Table->Delta_R_0;

  assert_right_model_definition( Table );

  R_W = b_W * e_Q;                           /* Threshold for Ants Persistence */

  #include <Model_Variables_Code.Include.c>

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {
 
    /* Under these conditions, four Dynamic Regimes at Stationarity are possible */
    if ( R_W < 1.0 ){
      printf("No ant population can persisten");
    }
  }
    else {
      printf("R_0 is only possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
      printf("the dynamics occurs in one single patch or local population\n");
      printf("But, here, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g)\n", 
	         Table->Lambda_C_0, Table->Lambda_R_0);
      printf("or the total number of patches is larger than 1 (N = %d)\n", Table->No_of_CELLS);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
  }

  return(R_W);
}




