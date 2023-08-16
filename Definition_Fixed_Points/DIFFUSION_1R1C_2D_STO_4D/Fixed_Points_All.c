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
    
    if ( Coexistence_Condition( Table ) == 0 ) {
      D = (Table->Lambda_R_0/Table->Beta_R - q)*(Table->Lambda_R_0/Table->Beta_R - q) + 4.0*Table->Lambda_R_0/Table->Beta_R;
      
      if (q > 0.0 ) {
	Vector_Stationarity_Lower[R]   = (double)Table->K_R/2.0 * (q - Table->Lambda_R_0/Table->Beta_R + sqrt(D));
	Vector_Stationarity_Lower[A]   = 0.0;
	Vector_Stationarity_Lower[RA]  = 0.0; 
	Vector_Stationarity_Lower[ARA] = 0.0;
      }
      else {
	if( Table->Lambda_R_0 == 0.0 ) {
	  Vector_Stationarity_Lower[R]   = 0.0;
	  Vector_Stationarity_Lower[A]   = 0.0;
	  Vector_Stationarity_Lower[RA]  = 0.0; 
	  Vector_Stationarity_Lower[ARA] = 0.0;
	}
	else {
	  Vector_Stationarity_Lower[R]   = (double)Table->K_R/2.0 * (q - Table->Lambda_R_0/Table->Beta_R + sqrt(D));
	  Vector_Stationarity_Lower[A]   = 0.0;
	  Vector_Stationarity_Lower[RA]  = 0.0; 
	  Vector_Stationarity_Lower[ARA] = 0.0;
	}
      }
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
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }
} 

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool; 
  double q, q_Star;
  
  q_Star = Table->Delta_C_0/Table->Alpha_C_0 + Table->Lambda_R_0/Table->Beta_R * (1.0 - Table->Alpha_C_0/Table->Delta_C_0);
  
  q = 1.0 - Table->Delta_R_0/Table->Beta_R;
  
  if (q < 0) {
    Condition_Bool = 0; /* No Resources: Extinction of Resources or Immigration-driven Dynamics if Lambda_R_0 is positive */
  }
  else { 
    if (q > q_Star) Condition_Bool = 1; /* Coexistence (with self-maintained resources (q > 0) */
    else            Condition_Bool = 0; /* Only self-maintained Resources: 
					   (Extinction of Consumers) */
  }
  
  return(Condition_Bool); 
}

double Coexistence_Condition_Double ( Parameter_Table * Table )
{
  int Condition_Bool;
  double Condition_Double; 
  double q, q_Star;

  q_Star = Table->Delta_C_0/Table->Alpha_C_0 + Table->Lambda_R_0/Table->Beta_R * (1.0 - Table->Alpha_C_0/Table->Delta_C_0);

  q = 1.0 - Table->Delta_R_0/Table->Beta_R;

  if (q < 0) {
    Condition_Bool = 0; /* No Resources: Extinction of Resources or Immigration-driven Dynamics if Lambda_R_0 is positive */
  }
  else { 
    if (q > q_Star) Condition_Bool = 1; /* Coexistence (with self-maintained resources (q > 0) */
    else            Condition_Bool = 0; /* Only self-maintained Resources: 
					   (Extinction of Consumers) */
  }
  
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
  
  if (Type_of_Stability == 0 ) {

    printf("%d: Stability: Only Resources\n", Type_of_Stability);
    
    Type_of_Stability_Double = (double)Type_of_Stability; 
    return(Type_of_Stability_Double);

  }
  else {
    
    double * Y0 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
    double * Y1 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
    double * Y2 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) ); 

    Fixed_Points_All( Table,
		      Y0, Y1, Y2, 1.0E-06);
    
    /* Eigen Values: Y1[k] + i Y2[k] from k=0 to Table->MODEL_STATE_VARIABLES-1 */
    GSL_Eigenvalue_Calculation(Y0, Table->MODEL_STATE_VARIABLES, Table,  Y1, Y2);
    // NR_Eigenvalue_Calculation(Y0, Table->MODEL_STATE_VARIABLES, Table,  Y1, Y2);

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
}
