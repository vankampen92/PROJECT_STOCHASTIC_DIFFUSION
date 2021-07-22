#include <MODEL.h>

void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon)
{
  int i; 
  double q, q_Star, D; 
  
  assert_right_model_definition( Table ); 

  #include <Model_Variables_Code.Include.c>
  
  if (Table->Lambda_C_0 == 0.0 && Table->No_of_CELLS == 1) {  

    q = 1.0 - Table->Delta_R_0/Table->Beta_R;
    
    if ( Coexistence_Condition( Table ) == 0 && q > 0.0 ) {

      D = (Table->Lambda_R_0/Table->Beta_R - q)*(Table->Lambda_R_0/Table->Beta_R - q) + 4.0*Table->Lambda_R_0/Table->Beta_R; 
      
      Vector_Stationarity_Lower[R]   = (double)Table->K_R/2.0 * (q - Table->Lambda_R_0/Table->Beta_R + sqrt(D));
      Vector_Stationarity_Lower[A]   = 0.0;
      Vector_Stationarity_Lower[RA]  = 0.0; 
      Vector_Stationarity_Lower[ARA] = 0.0; 
    }
    else {

      q_Star = Table->Delta_C_0/Table->Alpha_C_0 + Table->Lambda_R_0/Table->Beta_R * (1.0 - Table->Alpha_C_0/Table->Delta_C_0);
      
      Vector_Stationarity_Lower[R]   = (double)Table->K_R * Table->Delta_C_0/Table->Alpha_C_0;
      Vector_Stationarity_Lower[A]   = (double)Table->K_R * Table->Beta_R/Table->Alpha_C_0 * (q - q_Star); 
      Vector_Stationarity_Lower[RA]  = (double)Table->K_R * Table->Beta_R/Table->Alpha_C_0 * Table->Delta_C_0/Table->Nu_C_0 * (q - q_Star);
      if(Table->Chi_C_0 > 0.0)
	Vector_Stationarity_Lower[ARA] = (double)Table->K_R * Table->Beta_R/Table->Alpha_C_0*Table->Beta_R/Table->Alpha_C_0 * Table->Delta_C_0/Table->Nu_C_0 * Table->Chi_C_0/Table->Eta_C_0 * (q - q_Star)*(q - q_Star);
      else
	Vector_Stationarity_Lower[ARA] = 0.0; 
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
    printf("But, Lambda_C_0 is not zero (Lambda_C_0 = %g) or the number of patches\n", 
	   Table->Lambda_C_0);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);  
    Press_Key(); 
  }
} 

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool; 
  double q, q_Star;

  q_Star = Table->Delta_C_0/Table->Alpha_C_0 + Table->Lambda_R_0/Table->Beta_R * (1.0 - Table->Alpha_C_0/Table->Delta_C_0);

  q = 1.0 - Table->Delta_R_0/Table->Beta_R;

  if (q > q_Star) Condition_Bool = 1; /* Coexistence */
  else            Condition_Bool = 0; /* Only Resources: Extinction of Consumers */

  return(Condition_Bool); 
}

