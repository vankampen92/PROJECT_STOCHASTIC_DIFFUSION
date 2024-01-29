#include <MODEL.h>

void Fixed_Points_All( Parameter_Table * Table,
		                   double * Vector_Stationarity_Lower,
		                   double * Vector_Stationarity_Inter,
		                   double * Vector_Stationarity_Upper,
		                   double Epsilon)
{
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
      /* 0: Self-maintainance of ants by themselves is not possible */ 
      Vector_Stationarity_Lower[W]   = 0.0;
      Vector_Stationarity_Lower[Q]   = 0.0;
      Vector_Stationarity_Lower[F]   = 0.0;
      Vector_Stationarity_Lower[WF]  = 0.0;

      printf("Self-maintainance of the ant population is not possible\n");
    }
    else {

      // R_F = n_WF * R_W / (a_W * n_WF + d_W * (1.0 + n_WF) * e_Q * a_W); 
      
      R_F = R_W / (1.0 + (e_Q /a_F)*(1.0 + 1.0/n_WF));

      if ( R_F < 1.0 ) {
        /* 1: Only self-maintained resources with extintion of consumers */
	      Vector_Stationarity_Lower[W]   = (b_W - 1.0/e_Q);
	      Vector_Stationarity_Lower[Q]   = (b_W - 1.0/e_Q)/b_W;
	      Vector_Stationarity_Lower[F]   = 0.0;
	      Vector_Stationarity_Lower[WF]  = 0.0;

	      printf("Self-maintainance of the ant population by themselves in the absence\n");
        printf("of cosumers is possible, but consumers go extinct  (R_F < 1.0) \n");

      }
      else  {
        /* 2: Coexistence of resources and consumers */
	      xW  = (1.0 + n_WF)/a_F/n_WF;
        xQ  = e_Q * xW /(1.0 + e_Q*xW);
        xF  = (b_W*xQ - xW)/a_W/xW; 
        xWF = 1.0/n_F*xF; 
        
        Vector_Stationarity_Lower[W]   = xW; 
	      Vector_Stationarity_Lower[Q]   = xQ;
	      Vector_Stationarity_Lower[F]   = xF;
	      Vector_Stationarity_Lower[WF]  = xWF;

        printf("Peacefull coexistence of resources and consumers\n");
      }
  
      Vector_Stationarity_Lower[W]   *= Table->K_R;
      Vector_Stationarity_Lower[Q]   *= Table->K_R;
      Vector_Stationarity_Lower[F]   *= Table->K_R;
      Vector_Stationarity_Lower[WF]  *= Table->K_R;

      for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
        Vector_Stationarity_Inter[i] =                  Vector_Stationarity_Lower[i];
        Vector_Stationarity_Upper[i] =                  Vector_Stationarity_Lower[i];
        Table->Vector_Model_Variables_Stationarity[i] = Vector_Stationarity_Lower[i];
      }
    }
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g)\n", 
	         Table->Lambda_C_0, Table->Lambda_R_0);
    printf("or the total number of patches is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }
}

int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool;

  double b_W, e_Q, n_WF, n_F, a_W, a_F, d_W, R_W, R_F;
                                               /* Delta_W is  Delta_R_0 */
  b_W = Table->Beta_R/Table->Delta_R_0;        /* Delta_Q is  Delta_R_1 */
  e_Q = Table->Eta_R/Table->Delta_R_1;         /* Delta_F is  Delta_C_0 */
  n_WF = Table->Nu_C_0/Table->Delta_C_1;       /* Delta_WF is  Delta_C_1 */
  n_F  = Table->Nu_C_0/Table->Delta_C_0; 
  a_W  = Table->Alpha_C_0/Table->Delta_R_0;
  a_F  = Table->Alpha_C_0/Table->Delta_C_0;
  d_W  = Table->Delta_C_1/Table->Delta_R_0;

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {
    
    R_W = b_W * e_Q;                           /* Thresholds for Ants Persistence */
    R_F = R_W / (1.0 + (e_Q /a_F)*(1.0 + 1.0/n_WF));
    
    if( R_W > 1.0 && R_F > 1.0 ) Condition_Bool = 1;
    else                         Condition_Bool = 0;
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g)\n", 
	         Table->Lambda_C_0, Table->Lambda_R_0);
    printf("or the total number of patches is larger than 1 (N = %d)\n", Table->No_of_CELLS);
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

  double b_W, e_Q, n_WF, n_F, a_W, a_F, d_W, R_W, R_F;
                                               /* Delta_W is  Delta_R_0 */
  b_W = Table->Beta_R/Table->Delta_R_0;        /* Delta_Q is  Delta_R_1 */
  e_Q = Table->Eta_R/Table->Delta_R_1;         /* Delta_F is  Delta_C_0 */
  n_WF = Table->Nu_C_0/Table->Delta_C_1;       /* Delta_WF is  Delta_C_1 */
  n_F  = Table->Nu_C_0/Table->Delta_C_0; 
  a_W  = Table->Alpha_C_0/Table->Delta_R_0;
  a_F  = Table->Alpha_C_0/Table->Delta_C_0;
  d_W  = Table->Delta_C_1/Table->Delta_R_0;

  if (Table->Lambda_C_0 == 0.0 && Table->Lambda_R_0 == 0.0 && Table->No_of_CELLS == 1) {

    R_W = b_W * e_Q;           /* Thresholds for Ants Persistence */
    R_F = R_W / (1.0 + (e_Q /a_F)*(1.0 + 1.0/n_WF));
    
    if( R_W > 1.0 && R_F > 1.0 ) Condition_Bool = 1;
    else                         Condition_Bool = 0;
  }
  else {
    printf("Here, fixed points are only analytically possible if Lambda_C_0 and Lambda_R_0 are both 0 and \n");
    printf("the dynamics occurs in one single patch or local population\n");
    printf("But, Lambda_C_0 or Lambda_R_0 are not zero (Lambda_C_0 = %g and Lambda_R_0 = %g)\n", 
          Table->Lambda_C_0, Table->Lambda_R_0);
    printf("or the total number of patches is larger than 1 (N = %d)\n", Table->No_of_CELLS);
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

  double b_W, e_Q, n_WF, n_F, a_W, a_F, d_W, R_W, R_F;
                                               /* Delta_W is  Delta_R_0 */
  b_W = Table->Beta_R/Table->Delta_R_0;        /* Delta_Q is  Delta_R_1 */
  e_Q = Table->Eta_R/Table->Delta_R_1;         /* Delta_F is  Delta_C_0 */
  n_WF = Table->Nu_C_0/Table->Delta_C_1;       /* Delta_WF is  Delta_C_1 */
  n_F  = Table->Nu_C_0/Table->Delta_C_0; 
  a_W  = Table->Alpha_C_0/Table->Delta_R_0;
  a_F  = Table->Alpha_C_0/Table->Delta_C_0;
  d_W  = Table->Delta_C_1/Table->Delta_R_0;

  R_W = b_W * e_Q;           /* Thresholds for Ants Persistence */
  R_F = R_W / (1.0 + (e_Q /a_F)*(1.0 + 1.0/n_WF));
   
  Coexistence = Coexistence_Condition ( Table );
  
  if (Coexistence == 0 ) {

    if ( R_W < 1.0 )        Type_of_Stability = 0;
    /* 0: Self-maintainance of resources by themselves 
       in the absence of cosumers is not possible 
    */
    else if( R_W >= 1.0 )   Type_of_Stability = 1;
    /* 1: Only self-maintained resources 
       with extinction of consumers 
    */
    else{
      printf("Something very wrong in function Function_to_Type_of_Stability(...)\n");
      printf("R_W = %g\t R_F = %g\n", R_W, R_F);
	    
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
    Value_D = Y1[Index_Value_D];
    Index_D = Index_Value_D;
    
    if(Value_D < 0.0) {
      if(Y2[Index_D] == 0.0)       Type_of_Stability = 2;
      else if (Y2[Index_D] != 0.0) Type_of_Stability = 3;
    }
    else {
                                   Type_of_Stability = 4;
                                   // assert( Y2[Index_D] != 0.0 );
                                   assert(R_F > 1.0 );
      
    /* int X_apx_Y( double X, double Y, double ACU)
       
       When are two double variables, X and Y, approximately
       equal within some level of accuracy, ACU ? 

      // if (gsl_isnan(Y2[Index_D]) != 1) assert( Y2[Index_D] != 0.0 );
    */
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
