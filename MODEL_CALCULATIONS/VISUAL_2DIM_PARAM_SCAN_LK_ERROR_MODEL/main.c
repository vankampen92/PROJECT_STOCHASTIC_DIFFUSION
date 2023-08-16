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
   are defined as input arguments. 

   Compilation:
   
   . ~$ make MODEL=[TYPE_of_MODEL] 
  
   Execution:
                                                       
   . ~$ ./DIFFUSION_1R1C_2D -y0 4 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -H9 10.0 -I0 16 -m0 2.0 -M0 15.0 -A0 0.01 -d0 100 -H10 2.0 -I1 17 -m1 1.0 -M1 5.0 -A1 0.01 -d1 100 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 1 -D2 0 -a0 0 -Fn 1 -F0 Pseudo_Empirical_Data.dat -Y0 99 -tn 99 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000 -H4 1.0 -G30 R

   . ~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 12 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -H9 5.0 -I0 16 -m0 2.0 -M0 15.0 -A0 0.01 -d0 100 -H10 1.0 -I1 17 -m1 0.5 -M1 5.0 -A1 0.01 -d1 100 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 1 -D2 0 -a0 0 -tn 99 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000 -H4 1.5 -G30 R -Fn 1 -F0 Pseudo_Empirical_Data.dat -Y0 99

   . ~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 12 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -H9 5.0 -I0 16 -m0 2.0 -M0 15.0 -A0 0.01 -d0 100 -H10 1.0 -I1 17 -m1 0.5 -M1 5.0 -A1 0.01 -d1 100 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 1 -D2 0 -a0 0 -tn 99 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000 -H4 1.5 -G30 R -Fn 0 // Zero data files to read: no reading from data files: -Fn 0 to generate pseudo data //

-G30 R // Position of scale color bar: R, right side / L, left side / B, bottom side / T, top side .  
*/

gsl_rng * r; /* Global generator defined in main.c */
void Creating_Standard_Data_Matrix_from_Model ( Parameter_Table * Table,
						int i, 
						double **  Empirical_Data_Matrix ); 

void Minimum_Parameter_2D_Scan(Parameter_Table * Table,
			       int No_of_POINTS_1, int Input_Parameter_1,
			       int No_of_POINTS_2, int Input_Parameter_2,
			       double * W_GRID,
			       double * Likelihood_Minimum,
			       double * x_Value, 
			       double * y_Value); 

float * customized_contour_levels_0( Parameter_CPGPLOT * C )
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

float * customized_contour_levels_1( Parameter_CPGPLOT * C,
				     double Min_Value)
{
    int i;

    /* Two contour levels */
    C->NC = 2;
    float * clevels = (float *)calloc( C->NC, sizeof(float) );
    clevels[0] = Min_Value + 2.0;
    clevels[1] = Min_Value + 100.0;
    
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
  double Initial_Value_0, Initial_Value_1;
  
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
  printf("\n"); //Print_Press_Key(1,0,".");
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

  if( No_of_FILES > 0)
    pF = strcat(OBSERVED_DATA_FILE, Name_of_FILE[0]);          
  else 
    pF = strcat(OBSERVED_DATA_FILE, "Pseudo_Data_File.dat");          /* Default Name */

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
  if (TYPE_of_TIME_DEPENDENCE == 0) {      // -t4 0 (no time dependence)
                                           // -t4 1 (time-dependent parameters)
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
    					    1, Time_Empirical_Vector );
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

  if( No_of_FILES == 0) M_O_D_E_L___S_T_O( &Table ); /* Generation of Psedo Data */ 

  /* B E G I N : Observed Data Control Initization (Reserving Memmory)          */
  Observed_Data * Data = (Observed_Data *)calloc(1, sizeof(Observed_Data));
  Observed_Data_Alloc( Data, SUB_OUTPUT_VARIABLES, I_Time);
  printf(" Observed_Data structure has been correctly allocated\n");
  /* B E G I N : Reserving memmory for Parameter Fitting Structure */
  Parameter_Fitting * F = (Parameter_Fitting*)calloc(1,sizeof(Parameter_Fitting));
  F->Data                  = Data;
  F->Space                 = Space;
  F->Table                 = &Table;
  F->Minimization          = 0;     
  F->Bounded_Parameter_Set = 1;
  F->Function              = GSL_Function_to_Minimize_Error_Model; // GSL_Function_to_Minimize;
  /* double Function_to_Minimize( Parameter_Table * Table ) 
     This will be the function that will be past to the scanning function
     void generic_Function_Parameter_2Dim_Scan_Improved(Parameter_Table * Table, ...)
     The function 'Function_to_Minimize' is just wrapper for scanning the function 
     that the pointer 'F->Function' points to (see Optimization Library in ./Library 
     directory). 
  */
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
  
  int No_of_COLS;
  if( No_of_FILES > 0) No_of_COLS = F_y_GRID[0];  /* -Y0 [VALUE] */
  else                 No_of_COLS = I_Time;       /* -tn [VALUE] */

  Initial_Value_0 = AssignStructValue_to_VectorEntry(Input_Parameter_1, &Table);
  Initial_Value_1 = AssignStructValue_to_VectorEntry(Input_Parameter_2, &Table);  
  
  for (i = 0; i < Time.Realizations; i++ ) {

    if ( No_of_FILES > 0 )
      Reading_Standard_Data_Matrix_from_File( OBSERVED_DATA_FILE,
					      Empirical_Data_Matrix,
					      &SUB_OUTPUT_VARIABLES,
					      No_of_COLS, 
					      0, Name_of_Rows,
					      1, Time.Time_Vector );
    else
      Creating_Standard_Data_Matrix_from_Model ( &Table, i, 
						 Empirical_Data_Matrix );
    
    Writing_Standard_Data_Matrix( Empirical_Data_Matrix,
				  SUB_OUTPUT_VARIABLES, I_Time,
				  1, Name_of_Rows,
				  0, Time.Time_Vector);
    printf("Row Empirical Data Representation:\n");
    // Print_Press_Key(1,0,".");
    C_P_G___S_U_B___P_L_O_T_T_I_N_G___C_U_S_T_O_M_I_Z_E_D___T_I_T_L_E (&Table,
								       I_Time,
								       Time.Time_Vector,
								       Empirical_Data_Matrix,
								       0);
    // Print_Press_Key(1,0,".");
    
    Observed_Data_Initialization( Data, SUB_OUTPUT_VARIABLES,
				  I_Time, Empirical_Data_Matrix,
				  "" );
    printf(" Observed_Data structure has been correctly initiated with an instance\n");
    printf(" of (real or pseudo) emprical data \n");
    /*     E N D : -------------------------------------------------------------- */

    /* Back to input argument values in Table */
    AssignVectorEntry_to_Structure(&Table, Input_Parameter_1, Initial_Value_0);
    AssignVectorEntry_to_Structure(&Table, Input_Parameter_2, Initial_Value_1);
    
    P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( &Table, City_Par_Values) ;
    
    /* B E G I N : Main Function Call -------------------------------------------------*/
    
    double * W_GRID = (double *)malloc( No_of_POINTS_1 * No_of_POINTS_2 * sizeof(double) );
    /* 
       int Status =  generic_Function_Parameter_2Dim_Scan(&Table, 
       No_of_POINTS_1, Input_Parameter_1,
       No_of_POINTS_2, Input_Parameter_2,
       Function_to_Minimize, 
       W_GRID, "Negative LogLikelihood");
    */
    int X_LINEAR, Y_LINEAR;
    X_LINEAR = 0; Y_LINEAR = 0; /* Both axes linear */
    // X_LINEAR = 0; Y_LINEAR = 1; /* X axis linear and Y axis logarithmic */
    // X_LINEAR = 1; Y_LINEAR = 0; /* X axis logarithmic and Y axis linear */
    // X_LINEAR = 1; Y_LINEAR = 1; /* Both axes logarithmic */
    int Status =  generic_Function_Parameter_2Dim_Scan_Improved(&Table, 
								No_of_POINTS_1, Input_Parameter_1,
								No_of_POINTS_2, Input_Parameter_2,
								Function_to_Minimize, 
								W_GRID, "Negative LogLikelihood",
								X_LINEAR, Y_LINEAR);
    /*   E N D : ----------------------------------------------------------------------*/
    
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
      /* Comment out the following line if you want to write 
	 a predetermined title by the color wredge */
      Table.CPG->W_label[0] = '\0';
      
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
	 the function customized_contour_levels_[VALUE](...);
      */

      double Likelihood_Minimum, x_Val, y_Val;
      Minimum_Parameter_2D_Scan(&Table,
				No_of_POINTS_1, Input_Parameter_1,
				No_of_POINTS_2, Input_Parameter_2,
				W_GRID,
				&Likelihood_Minimum, &x_Val, &y_Val);
      
      //Table.CPG->contour_level = customized_contour_levels_0 ( Table.CPG );
      Table.CPG->contour_level = customized_contour_levels_1 ( Table.CPG,
							       Likelihood_Minimum );
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
      
      float x_Value = Parameter_Model_into_Vector_Entry(Input_Parameter_1, City_Par_Values);
      float y_Value = Parameter_Model_into_Vector_Entry(Input_Parameter_2, City_Par_Values);
      
      printf("%s=%f\t", Table.Symbol_Parameters[Input_Parameter_1], x_Value);
      printf("%s=%f\n", Table.Symbol_Parameters[Input_Parameter_2], y_Value);
      
      cpgslw(3);  /* Line width changing to 3     */
      cpgsci(12); /* Color Index changing to 12   */
      cpgpt1(x_Value, y_Value, 23);  /* Symbol 23 */
      
      float * xs = (float *)calloc(2, sizeof(float) );
      float * ys = (float *)calloc(2, sizeof(float) );
      xs[0] = 0.5* x_Value;  xs[1] = x_Value; /* A 40 % reduction */ 
      ys[0] = y_Value;       ys[1] = y_Value;
      
      cpg_XY_same_arrow( 2, xs, ys, 4, 1, 4);
      // cpg_XY_same_arrow( N, xs, ys, CPG->color_Index, CPG->type_of_Line, CPG->type_of_Width );
      
      free(xs);
      free(ys); 
#endif
 
      printf("Optimal Negative logLikelihood: %g\n", Likelihood_Minimum); 
      printf("%s=%f\t", Table.Symbol_Parameters[Input_Parameter_1], x_Val);
      printf("%s=%f\n", Table.Symbol_Parameters[Input_Parameter_2], y_Val);
      
      free (W_GRID);
      Print_Press_Key(1,0,"."); 
  }
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

void Minimum_Parameter_2D_Scan(Parameter_Table * Table,
				     int No_of_POINTS_1, int Input_Parameter_1,
				     int No_of_POINTS_2, int Input_Parameter_2,
				     double * W_GRID,
				     double * Likelihood_Minimum,
				     double * x_Value, 
				     double * y_Value)
{
	int n, k, j, i;
	int k_MIN, j_MIN; 
	double Minimum_Value, Value, Value_0, Value_1;
	
	Parameter_Space * S = Table->S;
	/* BEGIN : Allocating memory for saving data to plot a bifurcation  * * * * * * */
	/*         diagram for each variable  * * * * * * * * * * * * * * * * * * * * * */  
	double      ** z_SOL  = (double **)malloc( No_of_POINTS_2 * sizeof(double *) );
	for( i = 0; i < No_of_POINTS_2; i++){
	  z_SOL[i] = (double *)malloc( No_of_POINTS_1 * sizeof(double) );
	}
	double * x_Data  = (double *)malloc(No_of_POINTS_1 * sizeof(double) ); 
	double * y_Data  = (double *)malloc(No_of_POINTS_2 * sizeof(double) ); 
	/*   END : Allocating memory for saving dynamical data * * * * * */

	Minimum_Value = W_GRID[0]; 
	n = 0; 
	for( k = 0; k < No_of_POINTS_2; k++ ) {
  
	  Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min );
	  Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX );
	  
	  Value = Value_0 + k * (Value_1 - Value_0)/(double)(No_of_POINTS_2 - 1);
	  y_Data[k]= Value;
	  
	  Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min );
	  Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX );
	  
	  for( j = 0; j < No_of_POINTS_1; j++ ){
	    
	    Value = Value_0 + j * (Value_1 - Value_0)/(double)(No_of_POINTS_1 - 1);
	    
	    x_Data[j] = Value;
	    
	    z_SOL[k][j]    = W_GRID[n++]; 

	    Minimum_Value = MIN(Minimum_Value, z_SOL[k][j]);

	    if(Minimum_Value == z_SOL[k][j]) {
	      k_MIN = k;
	      j_MIN = j; 
	    }
	  }
	}

	//#if defined VERBOSE
	printf("Optimal Negative logLikelihood: %g\n", Minimum_Value); 
	printf("%s=%f\t", Table->Symbol_Parameters[Input_Parameter_1], x_Data[j_MIN]);
	printf("%s=%f\n", Table->Symbol_Parameters[Input_Parameter_2], y_Data[k_MIN]);
	//#endif

	* Likelihood_Minimum = Minimum_Value; 
	* x_Value = x_Data[j_MIN];
	* y_Value = y_Data[k_MIN];

	void Minimum_Parameter_2D_Scan(Parameter_Table * Table,
				     int No_of_POINTS_1, int Input_Parameter_1,
				     int No_of_POINTS_2, int Input_Parameter_2,
				     double * W_GRID,
				     double * Likelihood_Minimum,
				     double * x_Value, 
				     double * y_Value)
{
	int n, k, j, i;
	int k_MIN, j_MIN; 
	double Minimum_Value, Value, Value_0, Value_1;
	
	Parameter_Space * S = Table->S;
	/* BEGIN : Allocating memory for saving data to plot a bifurcation  * * * * * * */
	/*         diagram for each variable  * * * * * * * * * * * * * * * * * * * * * */  
	double      ** z_SOL  = (double **)malloc( No_of_POINTS_2 * sizeof(double *) );
	for( i = 0; i < No_of_POINTS_2; i++){
	  z_SOL[i] = (double *)malloc( No_of_POINTS_1 * sizeof(double) );
	}
	double * x_Data  = (double *)malloc(No_of_POINTS_1 * sizeof(double) ); 
	double * y_Data  = (double *)malloc(No_of_POINTS_2 * sizeof(double) ); 
	/*   END : Allocating memory for saving dynamical data * * * * * * * * * * * *  */

	Minimum_Value = W_GRID[0]; 
	n = 0; 
	for( k = 0; k < No_of_POINTS_2; k++ ) {
  
	  Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min );
	  Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX );
	  
	  Value = Value_0 + k * (Value_1 - Value_0)/(double)(No_of_POINTS_2 - 1);
	  y_Data[k]= Value;
	  
	  Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min );
	  Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX );
	  
	  for( j = 0; j < No_of_POINTS_1; j++ ){
	    
	    Value = Value_0 + j * (Value_1 - Value_0)/(double)(No_of_POINTS_1 - 1);
	    
	    x_Data[j] = Value;
	    
	    z_SOL[k][j]    = W_GRID[n++]; 

	    Minimum_Value = MIN(Minimum_Value, z_SOL[k][j]);

	    if(Minimum_Value == z_SOL[k][j]) {
	      k_MIN = k;
	      j_MIN = j; 
	    }
	  }
	}

	//#if defined VERBOSE
	printf("Optimal Negative logLikelihood: %g\n", Minimum_Value); 
	printf("%s=%f\t", Table->Symbol_Parameters[Input_Parameter_1], x_Data[j_MIN]);
	printf("%s=%f\n", Table->Symbol_Parameters[Input_Parameter_2], y_Data[k_MIN]);
	//#endif

	* Likelihood_Minimum = Minimum_Value; 
	* x_Value = x_Data[j_MIN];
	* y_Value = y_Data[k_MIN];
	
	/* BEGIN : Free memory!!!                                     * * * * * * */
	for( i = 0; i < No_of_POINTS_2; i++) free(z_SOL[i]); 
	free(z_SOL);
	free(x_Data); free(y_Data); 
	/*   END : Free memmory!!!                                    * * * * * * */
}
      
void Creating_Standard_Data_Matrix_from_Model ( Parameter_Table * Table,
						int i, 
						double **  Empirical_Data_Matrix )
{
  /* This is just to save the output variables corresponding to the i-th realization 
     generated by the function  M_O_D_EL___S_T_O( &Table ) into a straightforward 
     Empirical Data Matrix with the 'Psedo Data' that has been generated in the 
     i-th realization. 
  */
  int k,j;
  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++)
    for(j=0; j < Table->T->I_Time; j++)
	  Empirical_Data_Matrix[k][j] = Table->T->Variable[i][k][j];
  
}
