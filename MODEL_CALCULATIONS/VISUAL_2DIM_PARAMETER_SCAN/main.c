/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2021 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#include "global.h"

/* This code scans a function depending on model parameters as specified
   in Parameter_Table Structure. Any parametric function should be generically 
   defined as:
   
   double Function (Parameter_Table * Table); 

   For instance, a negative logLikelihood function, a Chi^2 Function, etc  

   Parameter ranges are defined in Parameter_Space structure. Parameters to scan 
   are defined as input arguments. For instance,  

   -sP 2 -I0 1 -H1 0.005 -m0 0.001  -M0 0.01    -d0 100     // Parameter: p_XY
         -I1 0 -H0 100.0 -m1 50.0 -M1 200.0     -d1 100     // Parameter: Beta_Y

   Compilation:
   
   . ~$ make MODEL=[TYPE_of_MODEL] 
   
   Execution:
                                                       
   . ~$ ./DIFFUSION_1R1C_2D -y0 4 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -sT 1.0E-06 -sN 300 -sP 2 -H9 10.0 -I0 16 -m0 0.0 -M0 50.0 -A0 0.01 -d0 10 -H10 2.0 -I1 17 -m1 0.0 -M1 10.0 -A1 0.01 -d1 10 -iP 0 -en 0 -eV 100.0 -DP 0 -DC 0 -D0 0 -D1 1 -D2 0 -a0 0 -Fn 1 -F0 Pseudo_Empirical_Data.dat -Y0 99 -tn 99 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000 -H4 1.0  

   Old example: 
   ./X2W2SILD-YSILD -y0 1  -sP 2 -I0 9 -m0 0.00001 -M0 0.001 -d0 400 -I1 0 -m1 50.0 -M1 120.0 -d1 400 -G0 1 -G1 1 -G14 R\\d\\fs0\\fn\\u -n 1 -v0 0 -en 0 -iP 0 
*/

gsl_rng * r; /* Global generator defined in main.c */


float * customized_contour_levels( Parameter_CPGPLOT * C )
{
    int i;

    /* Four contour levels */
    C->NC = 4;
    float * clevels = (float *)calloc( C->NC, sizeof(float) );
    clevels[0] = 1.0;
    clevels[1] = 2.5;
    clevels[2] = 5.0;
    clevels[3] =10.0; 
    
    return(clevels);
}

int main(int argc, char **argv)
{
  int i, k, key;
  char * pF; 
  Parameter_Table Table;
  Time_Control Time;
  Time_Dependence_Control Time_Dependence; 
  double Value_0, Value_1; 
  
  P_ARG = &Table;

#include "default.c"
  
  /* Extra default values */
  int No_of_POINTS_1    = 400;
  int Input_Parameter_1 = 12; 
  int No_of_POINTS_2    = 400;
  int Input_Parameter_2 = 25;

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  
  #include "include.Output_Variables.default.aux.c"
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C(   &Table );
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( &Table, Index_Output_Variables );
  printf(" Parameter_Table structure has been correctly allocated and initiated\n");

  Parameter_Model * City_Par_Values = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (&Table, City_Par_Values);
  printf(" Parameter_Model structure 'City_Par_Values' has been correctly allocated and initiated\n");
  
  /* B E G I N : Reserving memmory for Parameter Space */
  Parameter_Space * Error_Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  assert( No_of_ERROR_PARAMETERS == 0 );
  Table.E_Space = Error_Space;
 
  Parameter_Space * Initial_Condition_Space = (Parameter_Space *)calloc(1,
									sizeof(Parameter_Space));
  assert( No_of_IC == 0 );
  Table.IC_Space = Initial_Condition_Space;
  
  #include <include.Parameter_Space.default.aux.c>
  assert(No_of_PARAMETERS < Table.TOTAL_No_of_MODEL_PARAMETERS);
  Parameter_Space * Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  Parameter_Space_Alloc( Space, No_of_PARAMETERS, d);
  Parameter_Space_Initialization( Space, No_of_PARAMETERS, TOLERANCE, MAX_No_of_ITERATIONS,
				  d, Index, Ranges, Acc);
  Table.S = Space;
  printf("Parameter_Space structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */

  No_of_POINTS_1    = Space->N[0];
  Input_Parameter_1 = Space->Parameter_Index[0]; 
  No_of_POINTS_2    = Space->N[1];
  Input_Parameter_2 = Space->Parameter_Index[1];
  
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
  printf("\n"); //Press_Key();
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
  /* B E G I N : Observed Data Control Initization                                    */
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
  Press_Key(); 
  /* B E G I N :   Reserving memmory for Observed Data and Fitting Structure */
  Observed_Data * Data = (Observed_Data *)calloc(1, sizeof(Observed_Data));
  Observed_Data_Alloc( Data, SUB_OUTPUT_VARIABLES, I_Time);
  Observed_Data_Initialization( Data, SUB_OUTPUT_VARIABLES,
				I_Time, Empirical_Data_Matrix,
				"" );
  printf(" Observed_Data structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */
  
  /* B E G I N : Reserving memmory for Parameter Fitting Structure */
  Parameter_Fitting * F = (Parameter_Fitting*)calloc(1,sizeof(Parameter_Fitting));
  F->Data                  = Data;
  F->Space                 = Space;
  F->Table                 = &Table;
  F->Minimization          = 0;     
  F->Bounded_Parameter_Set = 1;
  F->Function              = GSL_Function_to_Minimize;
#if defined VERBOSE
  F->Verbose               = 1;     // 1: Verbose                // 0: Non Verbose
#else
  F->Verbose               = 0;     // 1: Verbose                // 0: Non Verbose
#endif
  F->Iteration             = 0;
  
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
  Table.Fitting_Data = (void *)F;
  
  printf("Parameter_Fitting structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */
  
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( &Table, City_Par_Values) ;
      
  double * W_GRID = (double *)malloc( No_of_POINTS_1 * No_of_POINTS_2 * sizeof(double) );
  int Status =  generic_Function_Parameter_2Dim_Scan(&Table, 
						     No_of_POINTS_1, Input_Parameter_1,
						     No_of_POINTS_2, Input_Parameter_2,
						     Function_to_Minimize, 
						     W_GRID, "Negative LogLikelihood");
#if defined CPGPLOT_REPRESENTATION
      /* BEGIN : 2D GRID cpgplot representation */
      /*********************************************************************/
      Table.CPG->X_label   = Table.Symbol_Parameters_Greek_CPGPLOT[Input_Parameter_1]; 
      Table.CPG->Y_label   = Table.Symbol_Parameters_Greek_CPGPLOT[Input_Parameter_2]; 
      /*********************************************************************/
      // Boundary(Input_Parameter_1, &Value_0, &Value_1);
      Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_min );
      Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_MAX );
      
      Table.CPG->ORIGIN_X    = Value_0;
      Table.CPG->X_Dimension = (Value_1 - Value_0);
      
      // Boundary(Input_Parameter_2, &Value_0, &Value_1);
      Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_min );
      Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_MAX );
      
      Table.CPG->ORIGIN_Y = Value_0;
      Table.CPG->Y_Dimension = (Value_1 - Value_0);
      
      Table.CPG->x_GRID  = No_of_POINTS_1; 
      Table.CPG->y_GRID  = No_of_POINTS_2;
      
      int Output_Variable  = Table.OUTPUT_VARIABLE_INDEX[0];
      Table.CPG->W_label   = Table.Output_Variable_Name[Output_Variable];

      Table.CPG->Title[0]='\0';  
      pF = strcat(Table.CPG->Title, "Negative Log Likelihood");
      
      int FIRST_PLOT = 0;
      double i_PLOT  = 0.0;
      C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___S_H_A_D_E_S( Table.CPG,
							    W_GRID, 
							    FIRST_PLOT,
							    Table.CPG->CPG_SCALE_W, 
							    Table.CPG->CPG_RANGE_W_0,
							    Table.CPG->CPG_RANGE_W_1,
							    i_PLOT );
      
      FIRST_PLOT = 1;
      Table.CPG->AUTOMATIC_CONTOUR = 0;
      /* If AUTOMATIC_CONTOUR is 0, the user should customized contours through
	 the function customized_contour_levels(...);
      */
      Table.CPG->contour_level = customized_contour_levels ( Table.CPG );
      C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___C_O_N_T_O_U_R( Table.CPG,
							      W_GRID, 
							      FIRST_PLOT,
							      Table.CPG->CPG_SCALE_W, 
							      Table.CPG->CPG_RANGE_W_0,
							      Table.CPG->CPG_RANGE_W_1,
							      i_PLOT );
      
      /* Annotating the countours by hand */
      // cpgptxt(float x, float y, float angle, float fjust,  const char *text);
      
      cpgslw(2);
      cpgptxt(0.00002, 62.0, 0.0, 0.0,  "1.0");
      cpgptxt(0.0001, 78.0, 0.0, 0.0,   "2.5");
      cpgptxt(0.0003, 92.0, 0.0, 0.0,  "5.0");
      
      /* Drawing arrow to emulate 0.6 reduction in 
	 p_YX transmission probability (from infecious female to male)
      */
      float x_Value = Parameter_Model_into_Vector_Entry(Input_Parameter_1, City_Par_Values);
      float y_Value = Parameter_Model_into_Vector_Entry(Input_Parameter_2, City_Par_Values);

      printf("%s=%f\t", Table.Symbol_Parameters[Input_Parameter_1], x_Value);
      printf("%s=%f\n", Table.Symbol_Parameters[Input_Parameter_2], y_Value);

      cpgslw(3);  /* Line width changing to 3     */
      cpgsci(12); /* Color Index changing to 12   */
      cpgpt1(x_Value, y_Value, 23);  /* Symbol 23 */
      
      float * xs = (float *)calloc(2, sizeof(float) );
      float * ys = (float *)calloc(2, sizeof(float) );
      xs[0] = x_Value;  xs[1] = 0.4 * x_Value; /* A 40 % reduction */ 
      ys[0] = y_Value;  ys[1] = y_Value;

      cpg_XY_same_arrow( 2, xs, ys, 4, 1, 4);
      // cpg_XY_same_arrow( N, xs, ys, CPG->color_Index, CPG->type_of_Line, CPG->type_of_Width );
      
      free(xs);
      free(ys); 
#endif  	
      
      free (W_GRID);

  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */ 
  Observed_Data_Free(Data); free(Data);

  if (TYPE_of_TIME_DEPENDENCE == 1)       // -t4 1
    Time_Dependence_Control_Free( &Time_Dependence, &Table );

#if defined CPGPLOT_REPRESENTATION
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif  

  #include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );
  free(Initial_Condition_Space);  
  free(Error_Space);
  free(F); 
 
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

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */
  
  printf("\nEnd of progam\n");
  return (0);
}


