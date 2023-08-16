/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2021 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#include "global.h"

gsl_rng * r; /* Global generator defined in main.c */

/* This code estimates ODE model parameters from data. Observed data are
   given as the result of a number of realizations of a stochastic model. 
   
   This main file will be then encapsulated as function of type:

   double Function_Relative_Accuracy( Parameter_Table * Table)

   which will give, as an output, a measure of the relative difference
   between an averge estimated value across realizations and the true value 
   (with which the stochastic realization has been obtained).

   This function function will allocate memory for parameter estimation
   (paramater space related estructures). This memmory will be freed 
   before leaving the function. 

   If there is time dependence (for instance, in just parameter 16: 
   (-t4 1 -DP 1 -DC- -D0 0 -D1 1 -D2 0 -P0 16), 
   then a time dependent file should be also provided through the input 
   agument list.

   Once the encapsulation works and the function: 
   
   double Function_Relative_Accuracy( Parameter_Table * Table)

   is ready to use, this function will be called repeatedly to explore
   the paremeter space in a main code: 

   / * B E G I N : Main Function Call -------------------------------------------------*/  
   // double * W_GRID = (double *)malloc( No_of_POINTS_1 * No_of_POINTS_2 * sizeof(double) );
   // int Status =  generic_Function_Parameter_2Dim_Scan(&Table, 
   //                                                    No_of_POINTS_1, Input_Parameter_1,
   // 						         No_of_POINTS_2, Input_Parameter_2,
   //  						         Function_Relative_Accuracy, 
   // 						         W_GRID, "Negative LogLikelihood");
   /*   E N D : -----------------------------------------------------------------------*/

   /* Compilation:

   . ~$ make

   Execution:

   . ~$ ./DIFFUSION_1R1C_2D -y0 4 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -sT 1.0E-08 -sN 300 -sP 2 -I0 16 -m0 2.0 -M0 15.0 -A0 0.01 -I1 17 -m1 1.0 -M1 5.0 -A1 0.01 -iP 0 -en 1 -e0 500.0 -em0 100 -eM0 1000  -DP 0 -DC 0 -D0 0 -D1 1 -D2 0 -a0 0 -Fn 1 -F0 Pseudo_Empirical_Data.dat -Y0 99 -tn 99 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000 -H4 1.0  
   
   (True parameters: -H9 [Alpha = 5.0] -H10 [Nu = 1.0] -H6 [Delta_C_0 = 1.0] used to generate pseudo data with DIFFUSION_1R1C) (Results file renamed as Full_Parameter_Set_Ordered.dat)

   According to the parameter correspondence in the Entropy paper, this would correspond to the following effective true parameters for the 2D system: 

   -H6 [Delta_C_0] approx 1  -H4 [Beta_R] approx 1  -H9 [Alpha] 10.0  -H10 [Nu] 2.0

`  . Pseudo_Data_File.dat has been generated with (in ./TEMPORAL_EVOLUTION)

   . ~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 4 -v0 0 -v1 1 -v2 2 -v3 12 -G0 2 -G1 2 -tn 100 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000

   Two important Notices: 
 1. -tn 100 and -Y0 100 should match!!! By contrast, -Y1 25 may or may not coincide, since 
    the funcion Time_Dependence_Control_Upload() uses splines to interpolate as many points 
    as necessary based on the points given in the file of time dependent parameters 
    (here Time_Dependent_Downloading_File.dat)
 2. -xn 1 (Random Intial Condition), then -xR 1 and -xN (total initial population) has to be
    non-zero. Otherwise, please, -xn 0 -xR 0  
*/

int main(int argc, char **argv)
{
  int i, j, k, n, s, z, key;
  double value, Min_Value, Data_Value, Theory_Value;
  double Likelihood_Value, Average_Likelihood_Value, Standard_Error_Value, Ave, Var;
  Parameter_Table Table;
  Time_Control Time;
  Time_Dependence_Control Time_Dependence;
  P_ARG = &Table;

#include "default.c"

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc, argv);

/* This is important to update Output Variable Index vector and dependent parameter 
   vector according to imput parameters from the command line 
*/
#include <include.Output_Variables.default.aux.c>

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C(   &Table );
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( &Table, Index_Output_Variables );
  printf(" Parameter_Table structure has been correctly allocated and initiated\n");

  Parameter_Model * Initial_Guess = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (&Table, Initial_Guess);
  printf(" Parameter_Model structure 'Initial_Guess' has been correctly allocated and initiated\n");
  Table.P = Initial_Guess; 
  
  /* B E G I N : Reserving memmory for Error Model Parameter Space */
  Parameter_Space * Error_Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  if( No_of_ERROR_PARAMETERS > 0 ) {
#include <include.Error_Control.default.aux.c>
    Parameter_Space_Alloc( Error_Space, No_of_ERROR_PARAMETERS, d_Error);
    Parameter_Table_Copy_into_Parameter_Model(Error_Space->Parameter_min, &Table);
    Model_Variables_Code_into_Parameter_Model(Error_Space->Parameter_min);
    Parameter_Table_Copy_into_Parameter_Model(Error_Space->Parameter_MAX, &Table);
    Model_Variables_Code_into_Parameter_Model(Error_Space->Parameter_MAX);
    Parameter_Table_Copy_into_Parameter_Model(Error_Space->Parameter_Acc, &Table);
    Model_Variables_Code_into_Parameter_Model(Error_Space->Parameter_Acc);
    Parameter_Space_Error_Initialization( Error_Space,
					  No_of_ERROR_PARAMETERS, d_Error,
					  Index_Error, Ranges_Error, Acc_Error );
    printf("Error Model Parameter_Space structure has been correctly allocated and initiated\n");
  }
  Table.E_Space = Error_Space;
  /*     E N D : ------------------------------------- */

  /* B E G I N : Reserving memmory for Initial Conditions Parameter Space */
  Parameter_Space * Initial_Condition_Space = (Parameter_Space *)calloc(1,
									sizeof(Parameter_Space));
  if( No_of_IC > 0 ) {
#include <include.Initial_Conditions.default.aux.c>
    Parameter_Space_Alloc( Initial_Condition_Space, No_of_IC, d_IC );
    Parameter_Table_Copy_into_Parameter_Model(Initial_Condition_Space->Parameter_min, &Table);
    Model_Variables_Code_into_Parameter_Model(Initial_Condition_Space->Parameter_min);
    Parameter_Table_Copy_into_Parameter_Model(Initial_Condition_Space->Parameter_MAX, &Table);
    Model_Variables_Code_into_Parameter_Model(Initial_Condition_Space->Parameter_MAX);
    Parameter_Table_Copy_into_Parameter_Model(Initial_Condition_Space->Parameter_Acc, &Table);
    Model_Variables_Code_into_Parameter_Model(Initial_Condition_Space->Parameter_Acc);
    Parameter_Space_IC_Initialization( Initial_Condition_Space,
				       No_of_IC, d_IC,
				       Index_IC, Ranges_IC, Acc_IC );
    printf("Initial Condition Parameter_Space structure has been correctly allocated and initiated\n");
  }
  Table.IC_Space = Initial_Condition_Space;
  /*     E N D : ------------------------------------- */

  /* B E G I N : Reserving memmory for Parameter Space */
#include <include.Parameter_Space.default.aux.c>
  if( No_of_PARAMETERS == Table.TOTAL_No_of_MODEL_PARAMETERS ) {
    /* Full parameter space is in place. See also Model_Variables_Code.c */
    for(i=0; i<Table.TOTAL_No_of_MODEL_PARAMETERS; i++) Index[i] = Table.Index[i];
    No_of_PARAMETERS = Table.TOTAL_No_of_MODEL_PARAMETERS;
  }
  Parameter_Space * Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  Parameter_Space_Alloc( Space, No_of_PARAMETERS, d);  // -sP 2
  Parameter_Space_Initialization( Space, No_of_PARAMETERS, TOLERANCE, MAX_No_of_ITERATIONS,
				  d, Index, Ranges, Acc);
  Table.S = Space;
  printf("Parameter_Space structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */

#if defined CPGPLOT_REPRESENTATION
  Table.CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  printf(" Parameter_CPGPLOT plotting structure has been correctly allocated and initiated\n");
#endif

#include <gsl_random_number_Setup.c>
#if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n");//Print_Press_Key(1,0,".");
  /*   END: Checking Random Number Generator Setup */
#endif

  char   ** Name_of_Rows          = (char **)calloc(SUB_OUTPUT_VARIABLES, sizeof(char *) );
  double ** Empirical_Data_Matrix = (double **)calloc( SUB_OUTPUT_VARIABLES, sizeof(double *) );
  for (i=0; i<SUB_OUTPUT_VARIABLES; i++) {
    key = Table.OUTPUT_VARIABLE_INDEX[i];
    Name_of_Rows[i]         = Table.Output_Variable_Symbol[key];
    Empirical_Data_Matrix[i] = (double *)calloc( I_Time, sizeof(double) );
  }

  /* B E G I N : Time Dependent Parameters, Observed Data, and Output files           */
  char * pF;
  char * TIME_PARAMETERS_FILE = (char *)calloc(1000, sizeof(char) ); /* Input files   */
  char * OBSERVED_DATA_FILE   = (char *)calloc(1000, sizeof(char) ); /* Input files   */
  char * PARAMETER_SET_FILE   = (char *)calloc(1000, sizeof(char) ); /* Output files  */

  FILE * DEMO;
  PARAMETER_SET_FILE[0] = '\0';
  pF = strcat(PARAMETER_SET_FILE, "Full_Parameter_Set.dat");
  DEMO = fopen(PARAMETER_SET_FILE, "w");

  /* In order for these default names to work properly, you need -Fn 0 in command line */
  OBSERVED_DATA_FILE[0] = '\0';
  pF = strcat(OBSERVED_DATA_FILE, "Observed_Data_File.dat");          /* Default Name */
  if( No_of_FILES > 0) strcpy(OBSERVED_DATA_FILE, Name_of_FILE[0]);   // -Fn 2

  TIME_PARAMETERS_FILE[0] = '\0';
  pF = strcat(TIME_PARAMETERS_FILE, "Time_Dependent_Parameters.dat"); /* Default Name */
  if( No_of_FILES > 1) strcpy(TIME_PARAMETERS_FILE, Name_of_FILE[1]);
  /*     E N D -----------------------------------------------------------------------*/

  /* B E G I N : Time Dependence Control common initization                           */
  double ** Matrix_of_COVARIATES; 
  double ** Type_1_Parameter_Values; 
  double  * Time_Empirical_Vector;
  char   ** Name_Rows_Dummy; 
  int No_of_EMPIRICAL_TIMES = F_y_GRID[1]; // No of Cols the time-dependent parameter file
  int No_of_Rows;                          // For example, -Y1 12 (see input argument list)
  if (TYPE_of_TIME_DEPENDENCE == 0) {      // -t4 1
    printf(" Time_Control structure will be allocated: \n");
    printf(" %d output variables of length %d points will be allocated\n",
	   SUB_OUTPUT_VARIABLES, I_Time);
    T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( &Time, &Table, I_Time);
    T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D( &Time, &Table, I_Time);
    printf(" Time_Control structure has been correctly allocated and set up\n");
  }
  else {
#include <include.Time_Dependence_Control.default.aux.c>
    printf(" Time_Dependence_Control and Time_Control structures will be allocated: \n");
    printf(" %d output variables of length %d points will be allocated\n",
	   SUB_OUTPUT_VARIABLES, I_Time);
    Time_Dependence_Control_Alloc(&Time, &Time_Dependence, &Table,
				  I_Time, TIME_DEPENDENT_PARAMETERS, No_of_COVARIATES);
    printf(" Both Time_Control and Time_Dependence_Control structures have been\n");
    printf(" correctly allocated\n");

    Type_1_Parameter_Values = (double **)calloc( TYPE_1_PARAMETERS, sizeof(double *));
    Time_Empirical_Vector   = (double * )calloc( No_of_EMPIRICAL_TIMES, sizeof(double));
    for(i = 0; i<TYPE_1_PARAMETERS; i++)
      Type_1_Parameter_Values[i] = (double *)calloc( No_of_EMPIRICAL_TIMES, sizeof(double));

    Reading_Standard_Data_Matrix_from_File( TIME_PARAMETERS_FILE,
    					    Type_1_Parameter_Values, &No_of_Rows,
					    No_of_EMPIRICAL_TIMES,
    					    0, Name_Rows_Dummy,
    					    1, Time_Empirical_Vector);
    assert( No_of_Rows == TYPE_1_PARAMETERS);

    Time_Dependence_Control_Upload_Optimized (&Time, &Time_Dependence, &Table,
    					      I_Time, No_of_EMPIRICAL_TIMES,
    					      TIME_DEPENDENT_PARAMETERS,
    					      TYPE_of_TIME_DEPENDENCE,        // -t4 1
    					      TYPE_0_PARAMETERS,              // -D0 0
    					      TYPE_1_PARAMETERS,              // -D1 1
    					      TYPE_2_PARAMETERS,              // -D2 0
    					      No_of_COVARIATES,               // -DC 0
    					      dependent_parameter, forcing_pattern, // -P0 16
    					      Matrix_of_COVARIATES,
    					      Type_1_Parameter_Values,
    					      Time_Empirical_Vector);
  }
  /*     E N D -----------------------------------------------------------------------*/
  
  int No_of_COLS = F_y_GRID[0]; // No of Columns in Observed Data File
  Reading_Standard_Data_Matrix_from_File( OBSERVED_DATA_FILE,
					  Empirical_Data_Matrix,
					  &SUB_OUTPUT_VARIABLES,
					  No_of_COLS, 
					  0, Name_of_Rows,
					  1, Time.Time_Vector );

  Writing_Standard_Data_Matrix( Empirical_Data_Matrix,
				SUB_OUTPUT_VARIABLES, I_Time,
				1, Name_of_Rows,
				0, Time.Time_Vector);
  Print_Press_Key(1,0,"."); 
  /* B E G I N :   Reserving memmory for Observed Data and Fitting Structure */
  Observed_Data * Data = (Observed_Data *)calloc(1, sizeof(Observed_Data));
  Observed_Data_Alloc( Data, SUB_OUTPUT_VARIABLES, I_Time);
  Observed_Data_Initialization( Data, SUB_OUTPUT_VARIABLES,
				I_Time, Empirical_Data_Matrix,
				"C-R Model" );
  printf(" Observed_Data structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */

  /* B E G I N :   Defining first row of "Full_Parameter_Set.dat" file
   */
  // fprintf(DEMO, "ID\t");
  for(i=0; i<Table.TOTAL_No_of_MODEL_PARAMETERS; i++) {
    key = Table.Index[i];
    fprintf(DEMO, "%s\t", Table.Symbol_Parameters[key]);
  }
  for(i=0; i<No_of_IC; i++) {
    key = Initial_Condition_Space->Parameter_Index[i];
    fprintf(DEMO, "%s\t", Table.Model_Variable_Symbol[key]);
  }
  for(i=0; i<No_of_ERROR_PARAMETERS; i++) {
    key = Error_Space->Parameter_Index[i];
    if(key < OUTPUT_VARIABLES_GENUINE_MAXIMUM)
      fprintf(DEMO, "Error(%s)\t", Table.Output_Variable_Symbol[key]);
    else
      fprintf(DEMO, "Error_Parameter_%d\t", i);
  }
  fprintf(DEMO, "NegLogLike\n");
  /*     E N D : --------------------------------------------------------
   */

  /* B E G I N : Reserving memmory for Fitting Structure -----------------
   */
  Parameter_Fitting * F = (Parameter_Fitting*)calloc(1,sizeof(Parameter_Fitting));
  Parameter_Fitting_Alloc( F, Realizations, &Table ); 
  Parameter_Fitting_Initialization(F, Realizations, Data, &Table);
  
  F->Minimization          = 1;     // 1: Function Minimization  // 0: Function Evaluation
  F->Bounded_Parameter_Set = 1;
  // F->Function              = GSL_Function_to_Minimize;
  F->Function              = GSL_Function_to_Minimize_Error_Model;
#if defined VERBOSE
  F->Verbose               = 1;     // 1: Verbose                // 0: Non Verbose
#else
  F->Verbose               = 0;     // 1: Verbose                // 0: Non Verbose
#endif
  /*     E N D : ----------------------------------------------------------
   */
  
  s = 0;
  int s_Attemps   = 0;   /* This counter will count number of random seeds */
  int Total_Tries = Realizations;                                // -tR 1000
  for(z=0; z<Realizations; z++) {  

    Parameter_Model_Copy_into_Parameter_Table(&Table, Initial_Guess);
    Random_Initial_Guess_within_Boundaries_Table(&Table, Space);
    Random_Error_Model_within_Boundaries_Table (&Table, Error_Space);
						
    /* Time_Dependence_Control_Upload_Optimized (&Time, &Time_Dependence, &Table,      */
    /* 					      I_Time, No_of_EMPIRICAL_TIMES,           */
    /* 					      TIME_DEPENDENT_PARAMETERS,               */
    /* 					      TYPE_of_TIME_DEPENDENCE,        // -t4 1 */
    /* 					      TYPE_0_PARAMETERS,              // -D0 0 */
    /* 					      TYPE_1_PARAMETERS,              // -D1 1 */
    /* 					      TYPE_2_PARAMETERS,              // -D2 0 */
    /* 					      No_of_COVARIATES,               // -DC 0 */
    /* 					      dependent_parameter, forcing_pattern, // -P0 16 */
    /* 					      Matrix_of_COVARIATES,                    */
    /* 					      Type_1_Parameter_Values,                 */
    /* 					      Time_Empirical_Vector);                  */

    if (TIME_DEPENDENT_PARAMETERS == 1) assert( Time_Dependence.No_of_TIMES == I_Time );
    
    /* Initial conditions from empirical data at the initial time (-xn 0 ) */
    /* Initial Condition from Empirical Data File (stored now in Empirical_Data_Matrix) */
    /* The system has only one cell (no spatial structure) */
    assert(Table.No_of_CELLS     == 1);
    assert(Table.No_of_RESOURCES == 1);
    Table.INITIAL_TOTAL_POPULATION = Empirical_Data_Matrix[0][0]; 
      
    F->Iteration             = z;
    
    if(F->Verbose == 1) {
      if(No_of_IC > 0)
	Parameter_Space_Write_Min_Max_IC_Values (Initial_Condition_Space, &Table );
      printf("\n");
      if(No_of_ERROR_PARAMETERS > 0)
	Parameter_Space_Write_Min_Max_Error_Values (Error_Space, &Table );
      printf("\n");
      Parameter_Space_Write_Min_Max_Values (Space, &Table );
      printf("\n");
    }

    printf("... Simplex bounded optimization starts right at this point within the\n");
    printf("... parameter space boundaries given above\n");
    // getchar();

    /* B E G I N :  This line of code invoques the optimizaton process for this random seed */
    Min_Value = GSL_Minimization_Driver( F );
    /*     E N D :  ----------------------------------------------------------------------- */

    if (F->Bounded_Parameter_Set == 1) {
      printf(" Min Value: NLL=%g\t Best Estimates: ", Min_Value);
      // fprintf(DEMO, "%d:\t", s_Attemps);

      for(i=0; i<Table.TOTAL_No_of_MODEL_PARAMETERS; i++) {
	key = Table.Index[i];
	value = AssignStructValue_to_VectorEntry(key, &Table);
	gsl_vector_set(F->Solution[s], i, value); 
	
	printf(" %s  = %g  ", Table.Symbol_Parameters[key], value);
	fprintf(DEMO, "%g\t", value);
      }
      for(i=0; i<No_of_IC; i++) {
	key = Initial_Condition_Space->Parameter_Index[i];
	value = Model_Variable_Initial_Condition_into_Vector_Entry_Table( key, &Table );
	gsl_vector_set(F->Solution[s], i+Table.TOTAL_No_of_MODEL_PARAMETERS, value);
	
	printf("%s = %g  ", Table.Model_Variable_Symbol[key], value);
	fprintf(DEMO, "%g\t", value);
      }
      for(i=0; i<No_of_ERROR_PARAMETERS; i++) {
	key = Error_Space->Parameter_Index[i];
	value = Error_Model_into_Vector_Entry_Table( key, &Table );
	gsl_vector_set(F->Solution[s], i+Table.TOTAL_No_of_MODEL_PARAMETERS+No_of_IC, value);
	
	if(key < OUTPUT_VARIABLES_GENUINE_MAXIMUM)
	  printf("%s = %g  ", Table.Output_Variable_Symbol[key], value);
	else
	  printf("Error_Parameter_%d = %g  ", i, value);
	fprintf(DEMO, "%g\t", value);
      }
      fprintf(DEMO, "%g\n", Min_Value);
      gsl_vector_set(F->Solution_Fitness, s, Min_Value); 
      
      s++; // Counting Successful attemts!!! 
    }
    printf("\n");

    s_Attemps++;

    printf(" Parameter Set Attemp No: %d\t Remaining Attemps: %d\t Successful Attemps: %d\n",
	   s_Attemps, Total_Tries-s_Attemps, s );
  }
  fclose(DEMO);

  Parametric_Configurations_from_Fitting_Structure_into_File (F,
							      "Full_Parameter_Set_Ordered.dat",
							      1);

  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */
  Parameter_Fitting_Free(F); free(F);
  
  Observed_Data_Free(Data); free(Data);

  if (TYPE_of_TIME_DEPENDENCE == 1)       // -t4 1
    Time_Dependence_Control_Free( &Time_Dependence, &Table );

#if defined CPGPLOT_REPRESENTATION
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif

  #include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );

  if( No_of_IC > 0 ) {
  #include <include.Initial_Conditions.default.free.c>
    Parameter_Space_Free(Initial_Condition_Space, No_of_IC);
  }
  free(Initial_Condition_Space);

  if( No_of_ERROR_PARAMETERS > 0 ) {
  #include <include.Error_Control.default.free.c>
    Parameter_Space_Free(Error_Space, No_of_IC);
  }
  free(Error_Space);

  #include <include.Output_Variables.default.free.c>

  #include <include.Time_Dependence_Control.default.free.c>
  if (TYPE_of_TIME_DEPENDENCE == 1) {                            // -t4 1
    for(i = 0; i<TYPE_1_PARAMETERS; i++) free(Type_1_Parameter_Values[i]);
    free(Type_1_Parameter_Values);
    free(Time_Empirical_Vector);
  }
  
  free(Name_of_Rows);
  for (i=0; i<SUB_OUTPUT_VARIABLES; i++)  free(Empirical_Data_Matrix[i]);
  free(Empirical_Data_Matrix);

  free(TIME_PARAMETERS_FILE);
  free(OBSERVED_DATA_FILE);
  free(PARAMETER_SET_FILE);

  free(Initial_Guess);

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */

  printf("\nEnd of progam\n");
  return (0);
}
