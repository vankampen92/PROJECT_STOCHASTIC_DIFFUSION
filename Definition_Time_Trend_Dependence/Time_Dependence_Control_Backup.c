#include <MODEL.h>

void Time_Dependence_Control_Alloc ( Time_Control * Time,
				     Time_Dependence_Control * TDC,
				     Parameter_Table * Table,
				     int No_of_TIMES,
				     int TIME_DEPENDENT_PARAMETERS,
				     int No_of_COVARIATES )
{
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

  T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C(Time, Table, No_of_TIMES);
  Table->T = Time; 
}

void Time_Dependence_Control_Upload (Time_Control * Time,
				     Time_Dependence_Control * TDC,
				     Parameter_Table * Table,
				     int No_of_TIMES,
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
  int i,j,k;
  int No_of_Rows;
  char ** Name_Dummy;
  double t;
  int * Function_Parameters = (int *)calloc(3, sizeof(int) ); 
  
  assert( TIME_DEPENDENT_PARAMETERS = TYPE_0_PARAMETERS + TYPE_1_PARAMETERS + TYPE_2_PARAMETERS);

  TDC->No_of_TIMES                  = No_of_TIMES;
  TDC->TIME_DEPENDENT_PARAMETERS    = TIME_DEPENDENT_PARAMETERS;
  TDC->TYPE_of_TIME_DEPENDENCE      = TYPE_of_TIME_DEPENDENCE;
  TDC->TYPE_0_PARAMETERS            = TYPE_0_PARAMETERS;
  TDC->TYPE_1_PARAMETERS            = TYPE_1_PARAMETERS;
  TDC->TYPE_2_PARAMETERS            = TYPE_2_PARAMETERS;
  TDC->No_of_COVARIATES             = No_of_COVARIATES;
  
  // TYPE_0_PARAMETERS is the number of parameters that, explicitly, depend
  // on COVARIATES. Covariates are time-varying variables that change at time a
  // increases. They may influence directly model parameters. We say they force
  // (externally) the system. This influence is usually mediated by response functions
  // that map these covariates into model parameters. 
  // TYPE_1_PARAMETERS are those time-dependent parametes that can be directly read from a file.
  // TYPE_2_PARAMETERS are those that functionally change with time according to certain given
  // functions. The type of functional dependency is determined by 'forcing_pattern'.  This
  // forcing pattern can be sinusoidal, sigmoidal, a linear increase, etc. 
  
  for(i = 0; i<TIME_DEPENDENT_PARAMETERS; i++) {
    TDC->Index_Dependent_Parameters[i] = dependent_parameter[i];
    TDC->Forcing_Pattern_Parameters[i] = forcing_pattern[i];
  }

  T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D(Time, Table, TDC->No_of_TIMES);

  if( TYPE_0_PARAMETERS > 0 ) {
    Reading_Standard_Data_Matrix_from_File( FILE_of_COVARIATES,
					    TDC->COVARIATES, &No_of_Rows, No_of_TIMES,
					    1, TDC->Name_of_COVARIATES,
					    1, TDC->Time_Vector);

    assert(No_of_Rows == No_of_COVARIATES);
    /* COVARIATES values determine Dependent_Parameters through response functions
     */
    // Upload_Dependent_Parameter_Values_from_Covariates (TDC, Table);
  }

  if(TYPE_1_PARAMETERS > 0 ) {
    double ** Type_1_Parameter_Values = (double **)calloc( TYPE_1_PARAMETERS, sizeof(double *));
    for(i = 0; i<TYPE_1_PARAMETERS; i++)
      Type_1_Parameter_Values[i] = (double *)calloc( TDC->No_of_TIMES, sizeof(double));

    Reading_Standard_Data_Matrix_from_File( FILE_of_TIME_PARAMETERS,
    					    Type_1_Parameter_Values, &No_of_Rows, No_of_TIMES,
    					    0, Name_Dummy,
    					    1, TDC->Time_Vector);

    assert( No_of_Rows == TYPE_1_PARAMETERS);

    for(i = TYPE_0_PARAMETERS; i < TYPE_0_PARAMETERS + TYPE_1_PARAMETERS; i++)
      for(j = 0; j < No_of_TIMES; j++)
	TDC->Dependent_Parameter[i][j] =  Type_1_Parameter_Values[i-TYPE_0_PARAMETERS][j];

    for(i = 0; i<TYPE_1_PARAMETERS; i++) free(Type_1_Parameter_Values[i]);
    free(Type_1_Parameter_Values);
  }

  if( TYPE_0_PARAMETERS > 0 || TYPE_1_PARAMETERS > 0 ) {
    Time->Time_0  = TDC->Time_Vector[0];
    Time->Time_1  = TDC->Time_Vector[No_of_TIMES-1];
    for(i = 0; i < TDC->No_of_TIMES; i++) Time->Time_Vector[i] = TDC->Time_Vector[i];
  }

  if(TYPE_2_PARAMETERS > 0) {
    for(i = 0; i < TYPE_2_PARAMETERS; i++) {
      k = i+TYPE_0_PARAMETERS+TYPE_1_PARAMETERS;
      for(j = 0; j<No_of_TIMES; j++) {
	t = Time->Time_Vector[j];
	TDC->Dependent_Parameter[k][j] =Time_Dependence_Resolve(Table,
								dependent_parameter[k],
								forcing_pattern[k], t);
      }
    }
  }

  Function_Parameters[0] = 23;
  Function_Parameters[1] = 24;
  Function_Parameters[2] = 25;
  
  Table->x_Bool = isDisjoint(Function_Parameters, 3,
			     Table->S->Parameter_Index, Table->S->No_of_PARAMETERS);   

  free(Function_Parameters);
  
  Table->TDC = TDC;
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
  
  free(TDC->Time_Vector);

  T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( P->T, P );
}
