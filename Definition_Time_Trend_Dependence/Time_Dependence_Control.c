#include <MODEL.h>

void Time_Dependence_Control_Alloc ( Time_Control * Time,
				     Time_Dependence_Control * TDC,
				     Parameter_Table * Table,
				     int No_of_TIMES, 
				     int TIME_DEPENDENT_PARAMETERS,
				     int No_of_COVARIATES )
{
  /* No_of_TIMES should be the total number of time intervals in which the numerical
     integration of the ODE system will be carried out. 
  */
  int i;

  TDC->Dependent_Parameter = (double **)calloc( TIME_DEPENDENT_PARAMETERS, sizeof(double * ));

  TDC->TIME_DEPENDENT_PARAMETERS    = TIME_DEPENDENT_PARAMETERS;
  for(i = 0; i < TIME_DEPENDENT_PARAMETERS; i++)
    TDC->Dependent_Parameter[i] = (double *)calloc(No_of_TIMES, sizeof(double));

  TDC->Index_Dependent_Parameters = (int *)calloc( TIME_DEPENDENT_PARAMETERS, sizeof(int));
  TDC->Forcing_Pattern_Parameters = (int *)calloc( TIME_DEPENDENT_PARAMETERS, sizeof(int));

  TDC->No_of_COVARIATES             = No_of_COVARIATES;
  if( No_of_COVARIATES > 0 ){ 
    TDC->COVARIATES         = (double **)calloc( No_of_COVARIATES, sizeof(double * ));
    TDC->Name_of_COVARIATES = (char **)calloc( No_of_COVARIATES, sizeof(char * ));
    for(i = 0; i < No_of_COVARIATES; i++) {
      TDC->COVARIATES[i] = (double *)calloc(No_of_TIMES, sizeof(double));
      TDC->Name_of_COVARIATES[i] = (char *)calloc(50, sizeof(char));
    }
  }
  
  TDC->Time_Vector = (double *)calloc(No_of_TIMES, sizeof(double));
    
  TDC->TDC_Function_Params = (double **)calloc(MODEL_PARAMETERS_MAXIMUM, sizeof(double *) );
  for(i = 0; i < MODEL_PARAMETERS_MAXIMUM; i++) 
    TDC->TDC_Function_Params[i] = (double *)calloc(No_of_TDC_FUNC_AUX_PARAM_MAX, sizeof(double));

  T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C(Time, Table, No_of_TIMES);
  Table->T = Time; 
}

void Time_Dependence_Control_Free ( Time_Dependence_Control * TDC, Parameter_Table * P)
{
  int i;

  for(i = 0; i < TDC->TIME_DEPENDENT_PARAMETERS; i++) free (TDC->Dependent_Parameter[i]);
  free (TDC->Dependent_Parameter);

  free(TDC->Index_Dependent_Parameters);
  free(TDC->Forcing_Pattern_Parameters);

  if( TDC->No_of_COVARIATES > 0 ){ 
    for(i = 0; i < TDC->No_of_COVARIATES; i++) {
      free(TDC->COVARIATES[i]);
      free(TDC->Name_of_COVARIATES[i]);
    }
    free(TDC->COVARIATES);
    free(TDC->Name_of_COVARIATES);
  }
  
  for(i = 0; i < MODEL_PARAMETERS_MAXIMUM; i++) 
    free(TDC->TDC_Function_Params[i]); 
  free(TDC->TDC_Function_Params); 
  
  free(TDC->Time_Vector);

  T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( P->T, P );
}

void Time_Dependence_Control_Upload (Time_Control * Time,
				     Time_Dependence_Control * TDC,
				     Parameter_Table * Table,
				     int No_of_TIMES,
				     int No_of_EMPIRICAL_TIMES, 
				     int TIME_DEPENDENT_PARAMETERS,
				     int TYPE_of_TIME_DEPENDENCE,
				     int TYPE_0_PARAMETERS,
				     int TYPE_1_PARAMETERS,
				     int TYPE_2_PARAMETERS,
				     int No_of_COVARIATES,
				     int * dependent_parameter,
				     int * forcing_pattern,
				     char * FILE_of_COVARIATES,
				     char * FILE_of_TIME_PARAMETERS)
{
  /* No_of_EMPIRICAL_TIMES is the number of times for which 
     we have data on changing parameter values. */ 
  int No_of_Rows; 
  int i,j,k;
  char ** Name_Dummy;
  double t;

  Table->TDC = TDC; /* First of all, the main structure (of type Parameter_Table), simply  
		       named as "Table" stores the pointer pointing to the 
		       Time Dependence Control structure named simply "TDC" */ 
  
  assert( TIME_DEPENDENT_PARAMETERS = TYPE_0_PARAMETERS + TYPE_1_PARAMETERS + TYPE_2_PARAMETERS);

  assert( No_of_EMPIRICAL_TIMES <= No_of_TIMES); 
  
  TDC->No_of_TIMES                  = No_of_TIMES;
  TDC->No_of_EMPIRICAL_TIMES        = No_of_EMPIRICAL_TIMES;
  TDC->TIME_DEPENDENT_PARAMETERS    = TIME_DEPENDENT_PARAMETERS;
  TDC->TYPE_of_TIME_DEPENDENCE      = TYPE_of_TIME_DEPENDENCE;
  TDC->TYPE_0_PARAMETERS            = TYPE_0_PARAMETERS;
  TDC->TYPE_1_PARAMETERS            = TYPE_1_PARAMETERS;
  TDC->TYPE_2_PARAMETERS            = TYPE_2_PARAMETERS;
  TDC->No_of_COVARIATES             = No_of_COVARIATES;
  
  // . TYPE_0_PARAMETERS is the number of parameters that, explicitly, depend
  // on COVARIATES. Covariates are time-varying variables that change as time a
  // increases. They may influence directly model parameters. We say they force
  // (externally) the system. This influence is usually mediated by response functions
  // that map these covariates into model parameters. 
  //
  // . TYPE_1_PARAMETERS are those time-dependent parametes that can be directly read from a file.
  //
  // . TYPE_2_PARAMETERS are those that functionally change with time according
  // to certain given functions. This dependency is direct (as opposed to mediated
  // by external (environmental parameters) and response functions). The type of functional
  // dependency is determined by the 'forcing_pattern'.  This forcing pattern can be
  // sinusoidal, sigmoidal, a linear increase, etc. 
  
  for(i = 0; i<TIME_DEPENDENT_PARAMETERS; i++) {
    TDC->Index_Dependent_Parameters[i] = dependent_parameter[i];
    TDC->Forcing_Pattern_Parameters[i] = forcing_pattern[i];
  }

  T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D(Time, Table, TDC->No_of_TIMES);

  // By default, the time vector of Time_Control structure should coincide with
  // the time vector of the Time_Dependence_Control structure:
  for(j = 0; j < No_of_TIMES; j++) TDC->Time_Vector[j] = Time->Time_Vector[j];

  // At the end of this function these two vectors should match!!! 
  
  if( TYPE_0_PARAMETERS > 0 ) {
    Reading_Standard_Data_Matrix_from_File( FILE_of_COVARIATES,
					    TDC->COVARIATES, &No_of_Rows, No_of_EMPIRICAL_TIMES,
					    1, TDC->Name_of_COVARIATES,
					    1, TDC->Time_Vector);

    assert(No_of_Rows == No_of_COVARIATES);
    /* COVARIATES values determine Dependent_Parameters through response functions
     */
    // Upload_Dependent_Parameter_Values_from_Covariates (TDC, Table);

    // Under construction... 
  }

  if(TYPE_1_PARAMETERS > 0 ) {
    double ** Type_1_Parameter_Values = (double **)calloc( TYPE_1_PARAMETERS, sizeof(double *));
    double  * Time_Empirical_Vector   = (double * )calloc( No_of_EMPIRICAL_TIMES, sizeof(double) );
    for(i = 0; i<TYPE_1_PARAMETERS; i++)
      Type_1_Parameter_Values[i] = (double *)calloc( No_of_EMPIRICAL_TIMES, sizeof(double));

    Reading_Standard_Data_Matrix_from_File( FILE_of_TIME_PARAMETERS,
    					    Type_1_Parameter_Values, &No_of_Rows,
					    No_of_EMPIRICAL_TIMES,
    					    0, Name_Dummy,
    					    1, Time_Empirical_Vector);
     
    assert( No_of_Rows == TYPE_1_PARAMETERS);

    if( TDC->No_of_EMPIRICAL_TIMES < TDC->No_of_TIMES ) {
      // Interpolation is required for each time-dependent parameter
      for(i = TYPE_0_PARAMETERS; i < TYPE_0_PARAMETERS + TYPE_1_PARAMETERS; i++) { 
	gsl_interp_accel * acc = gsl_interp_accel_alloc ();
	gsl_spline * spline    = gsl_spline_alloc (gsl_interp_cspline, No_of_EMPIRICAL_TIMES);
	gsl_spline_init (spline,
			 Time_Empirical_Vector, Type_1_Parameter_Values[i-TYPE_0_PARAMETERS],
			 No_of_EMPIRICAL_TIMES );

	for( j=0; j < No_of_TIMES; j++) 
	  TDC->Dependent_Parameter[i][j] = gsl_spline_eval( spline, Time->Time_Vector[j], acc); 
	  
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
      }

      for( j=0; j < No_of_TIMES; j++)  TDC->Time_Vector[j] = Time->Time_Vector[j];
    }
    else {

      if ( TDC->No_of_EMPIRICAL_TIMES > TDC->No_of_TIMES  ) { 
	printf(" Empirical No of Times should coincide with number of times\n");  
        printf(" in the numerical integration. Please make -tn [No of EMPIRICAL TIMES]\n");
      }
      
      assert( TDC->No_of_EMPIRICAL_TIMES == TDC->No_of_TIMES );
      
      for(i = TYPE_0_PARAMETERS; i < TYPE_0_PARAMETERS + TYPE_1_PARAMETERS; i++)
	for(j = 0; j < No_of_EMPIRICAL_TIMES; j++)
	  TDC->Dependent_Parameter[i][j] = Type_1_Parameter_Values[i-TYPE_0_PARAMETERS][j];

      for(j = 0; j < No_of_TIMES; j++) TDC->Time_Vector[j] = Time_Empirical_Vector[j];
    }
    
    for(i = 0; i<TYPE_1_PARAMETERS; i++) free(Type_1_Parameter_Values[i]);
    free(Type_1_Parameter_Values);
    free(Time_Empirical_Vector); 
  }

  if( TYPE_0_PARAMETERS > 0 || TYPE_1_PARAMETERS > 0 ) {
    Time->Time_0  = TDC->Time_Vector[0];
    Time->Time_1  = TDC->Time_Vector[No_of_TIMES-1];
    for(i = 0; i < TDC->No_of_TIMES; i++) Time->Time_Vector[i] = TDC->Time_Vector[i];
  }

  
  if(TYPE_2_PARAMETERS > 0) {

    // printf(" Time dependent parameters through time dependent functions:\n");
    Name_Dummy    = (char **)calloc(1, sizeof(char *) );
    Name_Dummy[0] = (char *)calloc(100, sizeof(char) );
    
    for(i = 0; i < TYPE_2_PARAMETERS; i++) {
      k = i+TYPE_0_PARAMETERS+TYPE_1_PARAMETERS;

      /* See Trend_Control.c for the definition of this function */
      Upload_Auxiliary_Parameter_Values_into_Table( Table,
						    dependent_parameter[k],
						    forcing_pattern[k] ); 
      Name_Dummy[0][0]= '\0';
      Symbol_to_Model_Parameters(dependent_parameter[k], Name_Dummy[0], Table); 

      // printf("%s\t", Name_Dummy[0]); 
      for(j = 0; j<No_of_TIMES; j++) {
	t = Time->Time_Vector[j];
	TDC->Dependent_Parameter[k][j] =Time_Dependence_Resolve(Table,
								dependent_parameter[k],
								forcing_pattern[k], t);
	// printf("%g, ", TDC->Dependent_Parameter[k][j] ); 	
      }
      // printf("\n");
    }

    // Print_Press_Key(1,0,".");

    free(Name_Dummy[0]);
    free(Name_Dummy); 
  }

  /* This is to see if certain parameters that control the way a true model parameter 
     depends on time have been included in Table->S->Parameter_Index, this is, in the
     list of parameters to be searched by the optimization algorithm 
  */

  /* This is to see if certain parameters that control the way a true model parameter 
     depends on time have been included in Table->S->Parameter_Index, this is, in the
     list of parameters to be searched by the optimization algorithm 
  */
  Table->x_Bool = 1;
  for(i=0; i<Table->S->No_of_PARAMETERS; i++) {
    if( Table->S->Parameter_Index[i] > MODEL_PARAMETERS_MAXIMUM) { 
      Table->x_Bool = 0;
      break;
    }
  }

  /* int * Function_Parameters = (int *)calloc(3, sizeof(int) );                     */
  /* Function_Parameters[0] = 23;                                                    */
  /* Function_Parameters[1] = 24;                                                    */
  /* Function_Parameters[2] = 25;                                                    */
  /* Table->x_Bool = isDisjoint(Function_Parameters, 3,                              */
  /* 			     Table->S->Parameter_Index, Table->S->No_of_PARAMETERS); */
  /* free(Function_Parameters);                                                      */
  /* Table->x_Bool is true if Function_Parameters and Table->S->Parameter_Index are 
     disjoint sets                                                                   */  
}

void Time_Dependence_Control_Upload_Optimized (Time_Control * Time,
					       Time_Dependence_Control * TDC,
					       Parameter_Table * Table,
					       int No_of_TIMES,
					       int No_of_EMPIRICAL_TIMES, 
					       int TIME_DEPENDENT_PARAMETERS,
					       int TYPE_of_TIME_DEPENDENCE,
					       int TYPE_0_PARAMETERS,
					       int TYPE_1_PARAMETERS,
					       int TYPE_2_PARAMETERS,
					       int No_of_COVARIATES,
					       int * dependent_parameter,
					       int * forcing_pattern,
					       double ** Matrix_of_COVARIATES,
					       double ** Type_1_Parameter_Values,
					       double  * Time_Empirical_Vector)
{
  /* No_of_EMPIRICAL_TIMES is the number of times for which 
     we have data on changing parameter values. */  
  int i,j,k;
  char ** Name_Dummy;
  double t;
  
  assert( TIME_DEPENDENT_PARAMETERS = TYPE_0_PARAMETERS + TYPE_1_PARAMETERS + TYPE_2_PARAMETERS);

  assert( No_of_EMPIRICAL_TIMES <= No_of_TIMES); 
  
  TDC->No_of_TIMES                  = No_of_TIMES;
  TDC->No_of_EMPIRICAL_TIMES        = No_of_EMPIRICAL_TIMES;
  TDC->TIME_DEPENDENT_PARAMETERS    = TIME_DEPENDENT_PARAMETERS;
  TDC->TYPE_of_TIME_DEPENDENCE      = TYPE_of_TIME_DEPENDENCE;
  TDC->TYPE_0_PARAMETERS            = TYPE_0_PARAMETERS;
  TDC->TYPE_1_PARAMETERS            = TYPE_1_PARAMETERS;
  TDC->TYPE_2_PARAMETERS            = TYPE_2_PARAMETERS;
  TDC->No_of_COVARIATES             = No_of_COVARIATES;
  
  // . TYPE_0_PARAMETERS is the number of parameters that, explicitly, depend
  // on COVARIATES. Covariates are time-varying variables that change as time a
  // increases. They may influence directly model parameters. We say they force
  // (externally) the system. This influence is usually mediated by response functions
  // that map these covariates into model parameters. 
  //
  // . TYPE_1_PARAMETERS are those time-dependent parametes that can be directly read from a file.
  //
  // . TYPE_2_PARAMETERS are those that functionally change with time according
  // to certain given functions. This dependency is direct (as opposed to mediated
  // by external (environmental parameters) and response functions). The type of functional
  // dependency is determined by the 'forcing_pattern'.  This forcing pattern can be
  // sinusoidal, sigmoidal, a linear increase, etc. 
  
  for(i = 0; i<TIME_DEPENDENT_PARAMETERS; i++) {
    TDC->Index_Dependent_Parameters[i] = dependent_parameter[i];
    TDC->Forcing_Pattern_Parameters[i] = forcing_pattern[i];
  }

  T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D(Time, Table, TDC->No_of_TIMES);

  // By default, the time vector of Time_Control structure should coincide with
  // the time vector of the Time_Dependence_Control structure:
  for(j = 0; j < No_of_TIMES; j++) TDC->Time_Vector[j] = Time->Time_Vector[j];

  // At the end of this function these two vectors should match!!! 
  
  if( TYPE_0_PARAMETERS > 0 ) {
   
    /* COVARIATES values determine Dependent_Parameters through response functions
     */
    // Upload_Dependent_Parameter_Values_from_Covariates (TDC, Table);

    // Under construction... 
  }

  if(TYPE_1_PARAMETERS > 0 ) {
   
    if( TDC->No_of_EMPIRICAL_TIMES < TDC->No_of_TIMES ) {
      // Interpolation is required for each time-dependent parameter
      for(i = TYPE_0_PARAMETERS; i < TYPE_0_PARAMETERS + TYPE_1_PARAMETERS; i++) { 
	gsl_interp_accel * acc = gsl_interp_accel_alloc ();
	gsl_spline * spline    = gsl_spline_alloc (gsl_interp_cspline, No_of_EMPIRICAL_TIMES);
	gsl_spline_init (spline,
			 Time_Empirical_Vector, Type_1_Parameter_Values[i-TYPE_0_PARAMETERS],
			 No_of_EMPIRICAL_TIMES );

	for( j=0; j < No_of_TIMES; j++) 
	  TDC->Dependent_Parameter[i][j] = gsl_spline_eval( spline, Time->Time_Vector[j], acc); 
	  
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
      }

      for( j=0; j < No_of_TIMES; j++)  TDC->Time_Vector[j] = Time->Time_Vector[j];
    }
    else {

      if ( TDC->No_of_EMPIRICAL_TIMES > TDC->No_of_TIMES  ) { 
	printf(" Empirical No of Times should coincide with number of times\n");  
        printf(" in the numerical integration. Please make -tn [No of EMPIRICAL TIMES]\n");
      }
      
      assert( TDC->No_of_EMPIRICAL_TIMES == TDC->No_of_TIMES );
      
      for(i = TYPE_0_PARAMETERS; i < TYPE_0_PARAMETERS + TYPE_1_PARAMETERS; i++)
	for(j = 0; j < No_of_EMPIRICAL_TIMES; j++)
	  TDC->Dependent_Parameter[i][j] = Type_1_Parameter_Values[i-TYPE_0_PARAMETERS][j];

      for(j = 0; j < No_of_TIMES; j++) TDC->Time_Vector[j] = Time_Empirical_Vector[j];
    } 
  }

  if( TYPE_0_PARAMETERS > 0 || TYPE_1_PARAMETERS > 0 ) {
    Time->Time_0  = TDC->Time_Vector[0];
    Time->Time_1  = TDC->Time_Vector[No_of_TIMES-1];
    for(i = 0; i < TDC->No_of_TIMES; i++) Time->Time_Vector[i] = TDC->Time_Vector[i];
  }

  if(TYPE_2_PARAMETERS > 0) {
    for(i = 0; i < TYPE_2_PARAMETERS; i++) {
      k = i+TYPE_0_PARAMETERS+TYPE_1_PARAMETERS;

      /* See Trend_Control.c for the definition of this function */
      Upload_Auxiliary_Parameter_Values_into_Table( Table,
						    dependent_parameter[k],
						    forcing_pattern[k] );
      
      for(j = 0; j<No_of_TIMES; j++) {
	t = Time->Time_Vector[j];
	TDC->Dependent_Parameter[k][j] =Time_Dependence_Resolve(Table,
								dependent_parameter[k],
								forcing_pattern[k], t);
	
      }
    }
  }

  /* This is to see if certain parameters that control the way a true model parameter 
     depends on time have been included in Table->S->Parameter_Index, this is, in the
     list of parameters to be searched by the optimization algorithm 
  */
  Table->x_Bool = 1;
  for(i=0; i<Table->S->No_of_PARAMETERS; i++) {
    if( Table->S->Parameter_Index[i] > MODEL_PARAMETERS_MAXIMUM) { 
      Table->x_Bool = 0;
      break;
    }
  }

  /* int * Function_Parameters = (int *)calloc(3, sizeof(int) ); */
  /* Function_Parameters[0] = 23;                                */
  /* Function_Parameters[1] = 24;                                */
  /* Function_Parameters[2] = 25;                                */
  /* Table->x_Bool = isDisjoint(Function_Parameters, 3,          */
  /* 			     Table->S->Parameter_Index, Table->S->No_of_PARAMETERS); */
  /* free(Function_Parameters);                                                      */
  /* Table->x_Bool is true if Function_Parameters and Table->S->Parameter_Index are 
     disjoint sets */  

  Table->TDC = TDC;
}
