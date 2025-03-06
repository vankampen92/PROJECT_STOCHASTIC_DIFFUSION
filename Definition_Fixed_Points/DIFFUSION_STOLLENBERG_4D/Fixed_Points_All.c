#include <MODEL.h>

void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon)
{
  int i;
  double d_P, b_R, bR_C, b_C, d_C, a, u, e_R;
  double q, w, x, y, z, x_C;
  double B_R;

  d_P = Table->Beta_R/Table->Delta_R_0;
  b_R = Table->Beta_R/Table->Delta_R_0;
  b_C = Table->Beta_C/Table->Delta_R_0;
  d_C = Table->Delta_C_0/Table->Delta_R_0;
  a   = Table->Alpha_C_0/Table->Delta_R_0;
  u   = Table->Nu_C_0/Table->Delta_R_0;
  e_R = Table->Eta_R/Table->Delta_R_0;

  assert_right_model_definition( Table );

  #include <Model_Variables_Code.Include.c>

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->Lambda_R_1 == 0.0 && Table->No_of_CELLS == 1) {

    /* Under these conditions, four Dynamic Regimes at Stationarity are possible */
    if ( b_R < (1.0 + d_P/e_R) ){
      /* 0: Self-maintainance of resources by themselves in the absence of cosumers is not possible */

      Vector_Stationarity_Lower[RP]  = 0.0;
      Vector_Stationarity_Lower[R]   = 0.0;
      Vector_Stationarity_Lower[A]   = 0.0;
      Vector_Stationarity_Lower[RA]  = 0.0;

      printf("Self-maintainance of resources by themselves in the absence of cosumers is not possible\n");
    }
    else {

      x_C  = (d_C + u)/a *1.0/(b_C/d_C - 1.0);
      B_R  = 1.0 - d_P/(b_R-1.0)/e_R;
      bR_C = 1.0 + d_P/e_R * 1.0 /(1-x_C); 
      
      if( b_C/d_C < 1.0 ) {
	/* 1: Only self-maintained resources with extintion of consumers */

	Vector_Stationarity_Lower[RP]  = b_R * B_R/(d_P + (1.0-B_R)*e_R);
	Vector_Stationarity_Lower[R]   = B_R;
	Vector_Stationarity_Lower[A]   = 0.0;
	Vector_Stationarity_Lower[RA]  = 0.0;

	printf("Self-maintainance of resources by themselves in the absence of cosumers is possible,\n");
	printf("but consumers go extinct  (b_C/d_C < 1.0) \n");

      }
      else if (x_C > 1.0) {

	Vector_Stationarity_Lower[RP]  = b_R * B_R/(d_P + (1.0-B_R)*e_R);
	Vector_Stationarity_Lower[R]   = B_R;
	Vector_Stationarity_Lower[A]   = 0.0;
	Vector_Stationarity_Lower[RA]  = 0.0;

	printf("Self-maintainance of resources by themselves in the absence of cosumers is possible,\n");
	printf("but consumers go extinct (x_C > 1.0)\n");
      }
      else if (b_R < bR_C) {

	Vector_Stationarity_Lower[RP]  = b_R * B_R/(d_P + (1.0-B_R)*e_R);
	Vector_Stationarity_Lower[R]   = B_R;
	Vector_Stationarity_Lower[A]   = 0.0;
	Vector_Stationarity_Lower[RA]  = 0.0;

	printf("Self-maintainance of resources by themselves in the absence of cosumers is possible,\n");
	printf("but consumers go extinct (x_C > 1.0)\n");
      }
      else { /* Well adapted animals (x_C < 1.0) and b_C/d_C > 1.0 */

	/* 2: Coexistence of resources and consumers */
	Vector_Stationarity_Lower[RP]  = b_R * x_C/(d_P + (1.0-x_C)*e_R);
	Vector_Stationarity_Lower[R]   = x_C;
	Vector_Stationarity_Lower[A]   = (b_R-1.0)*e_R/a*  (B_R - x_C)/(d_P + (1.0-x_C)*e_R);
	Vector_Stationarity_Lower[RA]  = (b_R-1.0)*e_R/a*1.0/(b_C/d_C - 1.0)* (B_R - x_C)/(d_P + (1.0-x_C)*e_R);

	printf("Peacefull coexistence of resources and consumers\n");

      }
    }

    Vector_Stationarity_Lower[RP]  *= Table->K_R;
    Vector_Stationarity_Lower[R]   *= Table->K_R;
    Vector_Stationarity_Lower[A]   *= Table->K_R;
    Vector_Stationarity_Lower[RA]  *= Table->K_R;

    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  Vector_Stationarity_Lower[i];
      Vector_Stationarity_Upper[i] =                  Vector_Stationarity_Lower[i];
      Table->Vector_Model_Variables_Stationarity[i] = Vector_Stationarity_Lower[i];
    }
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 or Lambda_R_1 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g and Lambda_R_1 = %g) or the number of patches\n",
	   Table->Lambda_C_0, Table->Lambda_R_0, Table->Lambda_R_1);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }
}

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool;
  double d_P, b_R, b_C, d_C, a, u, e_R;
  double x_C, bR_C0, bR_C1;

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->Lambda_R_1 == 0.0 && Table->No_of_CELLS == 1) {
    d_P = Table->Beta_R/Table->Delta_R_0;
    b_R = Table->Beta_R/Table->Delta_R_0;
    b_C = Table->Beta_C/Table->Delta_R_0;
    d_C = Table->Delta_C_0/Table->Delta_R_0;
    a   = Table->Alpha_C_0/Table->Delta_R_0;
    u   = Table->Nu_C_0/Table->Delta_R_0;
    e_R = Table->Eta_R/Table->Delta_R_0;

    x_C = (d_C + u)/a *1.0/(b_C/d_C - 1.0);

    bR_C0 = 1.0 + d_P/e_R;

    bR_C1 = 1.0 + d_P/e_R * 1.0 /(1-x_C); 

    if( b_R > bR_C0 && b_C/d_C > 1.0 && x_C < 1.0 && b_R > bR_C1) Condition_Bool = 1;

    else                                                          Condition_Bool = 0;

  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0, Lambda_R_0, and Lambda_R_1 are all 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
		printf("But, Lambda_C_0 or Lambda_R_0 or Lambda_R_1 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g and Lambda_R_1 = %g) or the number of patches\n",
		 Table->Lambda_C_0, Table->Lambda_R_0, Table->Lambda_R_1);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }

  return(Condition_Bool);
}

double Coexistence_Condition_Double ( Parameter_Table * Table )
{
  int Condition_Bool;
  double Condition_Double;

  double d_P, b_R, b_C, d_C, a, u, e_R;
  double x_C, bR_C0, bR_C1;

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->Lambda_R_1 == 0.0 && Table->No_of_CELLS == 1) {

    d_P = Table->Beta_R/Table->Delta_R_0;
    b_R = Table->Beta_R/Table->Delta_R_0;
    b_C = Table->Beta_C/Table->Delta_R_0;
    d_C = Table->Delta_C_0/Table->Delta_R_0;
    a   = Table->Alpha_C_0/Table->Delta_R_0;
    u   = Table->Nu_C_0/Table->Delta_R_0;
    e_R = Table->Eta_R/Table->Delta_R_0;

    x_C = (d_C + u)/a *1.0/(b_C/d_C - 1.0);

    bR_C0 = 1.0 + d_P/e_R;

    bR_C1 = 1.0 + d_P/e_R * 1.0 /(1-x_C);
    
    if( b_R > bR_C0 && b_C/d_C > 1.0 && x_C < 1.0 && b_R > bR_C1 ) Condition_Bool = 1;

    else                                                           Condition_Bool = 0;

  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0, Lambda_R_0, and Lambda_R_1 are all 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
		printf("But, Lambda_C_0 or Lambda_R_0 or Lambda_R_1 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g and Lambda_R_1 = %g) or the number of patches\n",
		 Table->Lambda_C_0, Table->Lambda_R_0, Table->Lambda_R_1);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }

  Condition_Double = (double)Condition_Bool;

  return(Condition_Double);
}

double Calculate_Stability_Stationary_Point( Parameter_Table * Table )
{
 /*
  *  Input arguments:
  *
  *  . Table, a pointer to the main data structure controling all parameters
  *  of the excution.

  *  Output arguments:
  *
  *  . Type_of_Stability, the value that takes this function for parameters values
  *  as defined in input Table.
 */

  int i, k;
  int Type_of_Stability, Coexistence;
  int Index_Value_D, Index_Value_S;
  int Index_D;
  double Value_D;
  double Type_of_Stability_Double;

  double d_P, b_R, b_C, d_C, a, u, e_R;
  double x_C, bR_C0, bR_C1;  /* Thresholds */

  d_P = Table->Beta_R/Table->Delta_R_0;
  b_R = Table->Beta_R/Table->Delta_R_0;
  b_C = Table->Beta_C/Table->Delta_R_0;
  d_C = Table->Delta_C_0/Table->Delta_R_0;
  a   = Table->Alpha_C_0/Table->Delta_R_0;
  u   = Table->Nu_C_0/Table->Delta_R_0;
  e_R = Table->Eta_R/Table->Delta_R_0;

  Coexistence = Coexistence_Condition ( Table );

  x_C   = (d_C + u)/a *1.0/(b_C/d_C - 1.0);
  
  bR_C0 = 1.0 + d_P/e_R;
  
  bR_C1 = 1.0 + d_P/e_R * 1.0 /(1-x_C);  
  
  if (Coexistence == 0 ) {

    if ( b_R <= bR_C0 )        Type_of_Stability = 0;
    /* 0: Self-maintainance of resources by themselves 
       in the absence of cosumers is not possible 
    */
    else if( b_C/d_C <= 1.0 )  Type_of_Stability = 1;
    /* 1: Only self-maintained resources 
       with extinction of consumers 
    */
    else if( x_C >= 1.0 )      Type_of_Stability = 1;
    /* 1: Only self-maintained resources 
       with extinction of consumers 
    */
    else if( b_R <= bR_C1 )    Type_of_Stability = 1;
    /* 1: Only self-maintained resources 
       with extinction of consumers 
    */
    else{
      printf("Something very wrong in function Function_to_Type_of_Stability(...)\n");
      printf("x_C = %g\t bR = %g\t bR_C0 = %g\t bR_C1 = %g\t b_C = %g\t d_C = %g\n",
	     x_C, b_R, bR_C0, bR_C1, b_C, d_C);  
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
  }
  else {
    /* When there is coexistence, we distinguish 3 different situations:
       2: Coexistence to a fixed point with non-damped oscillations
       3: Coexistence to a fixed point with damped oscillations
       4: Coexistence to a an periodic atractor of limit cycles
    */
    
    double * Y0 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
    double * Y1 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
    double * Y2 = (double *)calloc(Table->MODEL_STATE_VARIABLES, sizeof(double) );
    
    Fixed_Points_All( Table, Y0, Y1, Y2, 1.0E-06);

    Stationary_Solution_Feasibility_Control ( Table ); 

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

    /* Looking only at the dominant... */
    // Value_D = Y1[Index_Value_D];
    // Index_D = Index_Value_D;
    
    if(Value_D < 0.0) {
      if(Y2[Index_D] == 0.0)       Type_of_Stability = 2;
      else if (Y2[Index_D] != 0.0) Type_of_Stability = 3;
    }
    else {
                                   Type_of_Stability = 4;
				   assert( Y2[Index_D] != 0.0 );
    }
    
#if defined VERBOSE
    if (Type_of_Stability == 2) printf("%d: Stability: Non-Damped Oscillations\n",
					 Type_of_Stability);
    if (Type_of_Stability == 3) printf("%d: Stability: Damped Oscillations\n",
				       Type_of_Stability);
    if (Type_of_Stability == 4) printf("%d: Unstability: Limits Cycles\n",
				       Type_of_Stability);
    
    Write_Parameter_Table( Table, Table->TOTAL_No_of_MODEL_PARAMETERS );
#endif

    free(Y0); free(Y1); free(Y2);
  }

  Type_of_Stability_Double = (double)Type_of_Stability;
  
  return(Type_of_Stability_Double);
}
