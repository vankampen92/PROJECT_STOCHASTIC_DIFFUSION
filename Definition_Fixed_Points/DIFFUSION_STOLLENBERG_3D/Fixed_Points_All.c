#include <MODEL.h>

void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon)
{
  int i; 
  double b_R, b_C, d_C, a, u;
  double q, x, y, z; 

  b_R = Table->Beta_R/Table->Delta_R_0; 
  b_C = Table->Beta_C/Table->Delta_R_0;
  d_C = Table->Delta_C_0/Table->Delta_R_0;
  a   = Table->Alpha_C_0/Table->Delta_R_0;
  u   = Table->Nu_C_0/Table->Delta_R_0; 
  
  assert_right_model_definition( Table );

  assert_positive_model_parameters (Table); 

  #include <Model_Variables_Code.Include.c>
  
  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {  
    
    /* Under these conditions, four Dynamic Regimes at Stationarity are possible */
    if ( b_R < 1.0 ) {
      /* 0: Self-maintainance of resources by themselves in the absence of cosumers is not possible */
      
      Vector_Stationarity_Lower[R]   = 0.0;
      Vector_Stationarity_Lower[A]   = 0.0;
      Vector_Stationarity_Lower[RA]  = 0.0;

      printf("Self-maintainance of resources by themselves in the absence of cosumers is not possible\n"); 
    }
    else {

      q = (b_R-1)/b_R; 
      x = (d_C + u)/a *1.0/(b_C/d_C - 1.0);
      
      if( b_C/d_C < 1.0 ) {
	/* 1: Only self-maintained resources with extintion of consumers */
      	
	Vector_Stationarity_Lower[R]   = 1.0 - 1.0/b_R; 
	Vector_Stationarity_Lower[A]   = 0.0; 
	Vector_Stationarity_Lower[RA]  = 0.0;

	printf("Only self-maintained resources with extintion of consumers\n");
      }
      else if (x > 1) {

	Vector_Stationarity_Lower[R]   = 1.0 - 1.0/b_R; 
	Vector_Stationarity_Lower[A]   = 0.0; 
	Vector_Stationarity_Lower[RA]  = 0.0;

	printf("Only self-maintained resources with extintion of consumers\n");
      }
      else if ( q < x) {

	Vector_Stationarity_Lower[R]   = 1.0 - 1.0/b_R; 
	Vector_Stationarity_Lower[A]   = 0.0; 
	Vector_Stationarity_Lower[RA]  = 0.0;

	printf("Only self-maintained resources with extintion of consumers\n");
      }
      else {
	/* Well adapted animals (x_C < 1.0) and b_C/d_C > 1.0 */
	/* 2: Coexistence of resources and consumers */
	 
	  y = b_R/a * ( q - x);
       
	  z = 1/(b_C/d_C - 1.0) * y; 

	  Vector_Stationarity_Lower[R]   = x; 
	  Vector_Stationarity_Lower[A]   = y; 
	  Vector_Stationarity_Lower[RA]  = z;

	  printf("Peacefull coexistence of resources and consumers\n"); 
      }
    }

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
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g) or the number of patches\n", 
	   Table->Lambda_C_0, Table->Lambda_R_0);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }
} 

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool; 
  double b_R, b_C, d_C, a, u;
  double q, x; 

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {  
    b_R = Table->Beta_R/Table->Delta_R_0; 
    b_C = Table->Beta_C/Table->Delta_R_0;
    d_C = Table->Delta_C_0/Table->Delta_R_0;
    a   = Table->Alpha_C_0/Table->Delta_R_0;
    u   = Table->Nu_C_0/Table->Delta_R_0; 
    
    q = (b_R-1.0)/b_R;
    
    x = (d_C + u)/a *1.0/(b_C/d_C - 1.0);

    if( b_R > 1.0 && b_C/d_C > 1.0 && q > x && x < 1 ) Condition_Bool = 1;
    
    else                                               Condition_Bool = 0;
  
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g) or the number of patches\n", 
	   Table->Lambda_C_0, Table->Lambda_R_0);
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

  double b_R, b_C, d_C, a, u;
  double q, x; 

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {  
    
    b_R = Table->Beta_R/Table->Delta_R_0; 
    b_C = Table->Beta_C/Table->Delta_R_0;
    d_C = Table->Delta_C_0/Table->Delta_R_0;
    a   = Table->Alpha_C_0/Table->Delta_R_0;
    u   = Table->Nu_C_0/Table->Delta_R_0; 
    
    q = (b_R-1.0)/b_R;

    x = (d_C + u)/a *1.0/(b_C/d_C - 1.0);

    if( b_R > 1.0 && b_C/d_C > 1.0 && q > x && x < 1) Condition_Bool = 1;
       
    else                                              Condition_Bool = 0;
   
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g) or the number of patches\n", 
	   Table->Lambda_C_0, Table->Lambda_R_0);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
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
   *  . Type_of_Stability, the value that takes this function for parameters values 
   *  as defined in input Table.  
*/
  
  int i, k;
  int Type_of_Stability, Coexistence;
  int Index_Value_D, Index_Value_S;
  int Index_D; 
  double Value_D;
  double Type_of_Stability_Double; 

  double b_R, b_C, d_C, a, u;
  double q, x, bR_C; 

  b_R = Table->Beta_R/Table->Delta_R_0; 
  b_C = Table->Beta_C/Table->Delta_R_0;
  d_C = Table->Delta_C_0/Table->Delta_R_0;
  a   = Table->Alpha_C_0/Table->Delta_R_0;
  u   = Table->Nu_C_0/Table->Delta_R_0; 
    
  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {  

    Coexistence = Coexistence_Condition ( Table );

    q = (b_R-1.0)/b_R;

    x = (d_C + u)/a *1.0/(b_C/d_C - 1.0);

    bR_C = 1.0/(1.0 - x);  /* b_R < bR_C is equivalent to q < x */
    
    if (Coexistence == 0 ) {
      /* When there is no coexistence, we distinguish 3 different situations:
	 0: Self-maintainance of resources by themserlves in the absence of cosumers is not possible
	 1: Only self-maintained resources with extintion of consumers 
	 3: Over-exploitation of resources by greedy consumers and system total collapse 
      */
      if ( b_R <= 1.0 )             Type_of_Stability = 0; 
      /* 0: Self-maintainance of resources by themselves in the absence of cosumers is not possible */
      else if( b_C/d_C <= 1.0 )     Type_of_Stability = 1; 
      /* 1: Only self-maintained resources with extintion of consumers */
      else if( x >= 1.0 )           Type_of_Stability = 1; 
      /* 1: Only self-maintained resources with extintion of consumers */
      else if( q <= x )             Type_of_Stability = 1; 
      /* 1: Over-exploitation of resources by greedy consumers. Consumers go extinct and the 
	    system bounces back upto carrying capacity. In the end, only self-maintained 
            resources persist in the system*/
	      
      else {
	printf("Something very wrong in function Function_to_Type_of_Stability(...)\n"); 
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
      
      Fixed_Points_All( Table,
			Y0, Y1, Y2, 1.0E-06);

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
      
      if(Value_D < 0.0) {
	if(Y2[Index_D] == 0.0)       Type_of_Stability = 2;
	else if (Y2[Index_D] != 0.0) Type_of_Stability = 3;
      }
      else                           Type_of_Stability = 4; 

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
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g) or the number of patches\n", 
	   Table->Lambda_C_0, Table->Lambda_R_0);
    printf("is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }
  
  Type_of_Stability_Double = (double)Type_of_Stability;
  
  return(Type_of_Stability_Double);
}


  
