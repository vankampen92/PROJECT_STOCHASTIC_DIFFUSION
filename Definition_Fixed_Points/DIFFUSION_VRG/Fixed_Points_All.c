#include <MODEL.h>

void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon)
{
  int i; 
  double q, q_Star, a_Star, K_R, G_Star, B_Star; 
  
  assert_right_model_definition( Table ); 

  #include <Model_Variables_Code.Include.c>

  assert(Table->No_of_CELLS == 1);

  K_R = (double)Table->K_R; 
  
  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 ) {  

    q = 1.0 - Table->Delta_R_0/Table->Beta_R;

    assert(q > 0.0); 
    
    if ( Coexistence_Condition( Table ) == 0 && q > 0.0 ) {

      Vector_Stationarity_Lower[V]   = q * K_R; 
      Vector_Stationarity_Lower[R]   = 0.0;
      Vector_Stationarity_Lower[G]   = 0.0; 

    }
    else {

      
      q_Star = 1.0 - Table->Delta_C_0/Table->Theta_C/Table->Alpha_C_0; 
      a_Star = Table->Beta_C*q_Star-Table->Lambda_C_1-Table->Delta_C_1;
      
      Vector_Stationarity_Lower[V] = Table->Delta_C_0/Table->Theta_C/Table->Alpha_C_0 * K_R;

      G_Star  = (a_Star + sqrt(a_Star*a_Star + 4.0*Table->Lambda_C_1*q_Star*Table->Beta_C))/2.0/Table->Beta_C * K_R; 

      B_Star  = Table->Beta_R + Table->p_1 * G_Star; 

      Vector_Stationarity_Lower[R] = (B_Star * (q_Star-G_Star/K_R) - Table->Delta_R_0)/Table->Alpha_C_0 * K_R;

      Vector_Stationarity_Lower[G] = G_Star;
       
    }

    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  Vector_Stationarity_Lower[i];
      Vector_Stationarity_Upper[i] =                  Vector_Stationarity_Lower[i];
      Table->Vector_Model_Variables_Stationarity[i] = Vector_Stationarity_Lower[i];
    }
   
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 is 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, either Lambda_C_0 is not zero (Lambda_C_0 = %g) or the number of patches\n", 
	   Table->Lambda_C_0);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("All vectors at stationarity are set to zero\n");  
    Print_Press_Key(1,0,"."); 
  }
} 

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool; 
  double q_1, q_2, q_3, B_Star, G_Star, K_R;

  K_R  = (double)Table->K_R;

  q_1  = 1.0 - Table->Delta_C_0/Table->Alpha_C_0/Table->Theta_C; 

  q_2  = Table->Delta_C_1 - Table->Beta_C * q_1; 

  G_Star = Table->Lambda_C_1*q_1/(Table->Delta_C_1 - Table->Beta_C*q_1) * K_R;
  B_Star = Table->Beta_R + Table->p_1 * G_Star; 

  q_3  = q_1 - Table->Delta_R_0/B_Star;
    
  if ( q_1 > 0.0 && q_2 > 0.0 && q_3 > 0.0 ) Condition_Bool = 1; /* (V*, R*, G*) */
  else                                       Condition_Bool = 0; /*              */

  return(Condition_Bool); 
}

