#include <MODEL.h>

void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon)
{
  int i; 
  double q, q_Star, D;
  double K_R, y_R, Theta, A_0;

  assert_right_model_definition( Table ); 

  #include <Model_Variables_Code.Include.c>

  K_R   = (double)Table->K_R;
  y_R   = Table->TOTAL_No_of_RESOURCES; /* Resource density (in number of resource units) */
  Theta = Table->Alpha_C_0 * y_R/K_R;
  A_0   = (double)Table->TOTAL_No_of_CONSUMERS;
  
  if (Table->Lambda_R_0 == 0.0 && Table->Lambda_C_0 == 0.0 && Table->No_of_CELLS == 1) {  

      Vector_Stationarity_Lower[A]   = Table->Nu_C_0 / (Table->Nu_C_0 + Theta) * A_0;
      Vector_Stationarity_Lower[RA]  = Theta         / (Table->Nu_C_0 + Theta) * A_0;
    

    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  Vector_Stationarity_Lower[i];
      Vector_Stationarity_Upper[i] =                  Vector_Stationarity_Lower[i];
      Table->Vector_Model_Variables_Stationarity[i] = Vector_Stationarity_Lower[i];
    }
   
  }
  else {
    printf("Here, fixed points are only analytically possible if\n");
    printf("Lambda_C_0 is 0 Lambda_R_0 is zero, and, in addition,\n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 is not zero (Lambda_C_0 = %g)\n", Table->Lambda_C_0);
    printf(" or  Lambda_R_0 is not zero (Lambda_R_0 = %g)\n", Table->Lambda_R_0);
    printf(" or  the number of patches is larger than 1 (N = %d)\n", Table->No_of_CELLS);

    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }
} 

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool; 
  
  Condition_Bool = 1;
  
  return(Condition_Bool); 
}

double Coexistence_Condition_Double ( Parameter_Table * Table )
{
  int Condition_Bool;
  double Condition_Double; 

  Condition_Bool = 1;
  
  Condition_Double = (double)Condition_Bool;
  
  return(Condition_Double); 
}

double Function_to_Type_of_Stability( Parameter_Table * Table )
{
  /* 	
   *  Input arguments:
   *
   *  . Table, a pointer to the main data structure controling all parameters
   *  of the excution.
  
   *  Output arguments:
   *
   *  . value, the value that takes this function for parameters values 
   *  as defined in input Table. 
   
   In this example, this is a three-valued function. Return values are

   . 0, for only resources or no resources at all. 

   . 1, for stable coexistce of resources and consumers

   . 2, for stable coexistce of resources and consumers and damped oscillations 
*/
  
  int i, k;
  int Type_of_Stability;
  int Index_Value_D, Index_Value_S;
  int Index_D; 
  double Value_D;
  double Type_of_Stability_Double; 

  Type_of_Stability = Coexistence_Condition ( Table ); 
  
    
  double * Y0 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
  double * Y1 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
  double * Y2 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) ); 
  
  Fixed_Points_All( Table, Y0, Y1, Y2, 1.0E-06);
    
  /* Eigen Values: Y1[k] + i Y2[k] from k=0 to Table->MODEL_STATE_VARIABLES-1 */
 
  GSL_Eigenvalue_Calculation(Y0, Table->MODEL_STATE_VARIABLES, Table,  Y1, Y2);
  
  Dominant_Eigenvalue_Calculation(Y1, Y2, Table->MODEL_STATE_VARIABLES,
				  &Index_Value_D, &Index_Value_S); 
  
  if(Y1[Index_Value_D] == 0.0) {
    Value_D = Y1[Index_Value_S];
    Index_D = Index_Value_S;
  }
  else {
    Value_D = Y1[Index_Value_D];
    Index_D = Index_Value_D;
  }
  
  if(Value_D < 0.0) {
    if(Y2[Index_D] == 0.0)       Type_of_Stability = 1;
    else if (Y2[Index_D] != 0.0) Type_of_Stability = 2;
  }
  else Type_of_Stability = 3; 
  
  if (Type_of_Stability == 1) printf("%d: Stability: Non-Damped Oscillations\n",
				     Type_of_Stability); 
  if (Type_of_Stability == 2) { printf("%d: Stability: Damped Oscillations\n",
				       Type_of_Stability);
    // Print_Press_Key(1,0,".");
  }
  if (Type_of_Stability == 3) { printf("%d: Unstability: Limits Cycles\n",
				       Type_of_Stability);
    // Print_Press_Key(1,0,".");
  }
  
  Write_Parameter_Table( Table, Table->TOTAL_No_of_MODEL_PARAMETERS ); 
  // Print_Press_Key(1,0,".");
  
  free(Y0); free(Y1); free(Y2); 
  
  Type_of_Stability_Double = (double)Type_of_Stability; 
  return(Type_of_Stability_Double);
}
