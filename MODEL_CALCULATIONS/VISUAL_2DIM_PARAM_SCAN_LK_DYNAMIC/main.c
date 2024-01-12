/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2021 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#include "global.h"

/* This code scans a function depending on model parameters as specified
   in Parameter_Table Structure. Any parametric function should be generically 
   defined as:
   
   double Function (Parameter_Table * Table); 

   In this example, several negative logLikelihood functions can be analized. 
   In the directory ./VISUAL_2DIM_PARAM_SCAN_LK_DYNAMICS, main functions use 
   time-dependent likelihoods for different models (as opposed to likelihood
   functions that are defined at stationarity, which are addressed in 
   the directory ./VISUAL_2DIM_PARAM_SCAN_LK_STATIONARITY). The definitions
   of the likelihood functions themselves can be found in the directory
   ./Library/Optimization_Library/ where the functions are implemented 
   as the static library libda_Optimization.a 
 
   These definitions are based on the (in general) multivariate probability 
   distributions associated to a jump process. The temporal evolution of this 
   object is considered. For some models, there is a close expression for the 
   time evolution of this probability (MODELl=DIFFUSION_HII_nD), but most 
   models rely on the numerical integration of the master equation.    

   Parameter ranges (in 2D) are defined in Parameter_Space structure. Parameters 
   to scan are defined as input arguments. 

   Compilation:
   
   . ~$ make MODEL=[TYPE_of_MODEL] 
  
   Execution:
   
   PsedoData (output variables: -n 2 -v0 0 -v1 1): (n_A,n_RA) at stationarity for a number of realizations (-tR 100); 

   MODEL: BEDDINGTON-DeANGELIS. Scanning 2D Parameter Space (Alpha_C_0, Nu_C_0):
   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 3 -G1 3 -sT 1.0E-06 -sN 300 -sP 2 -H9 2.5 -I0 16 -m0 0.5 -M0 5.0 -A0 0.01 -d0 200 -H10 10.0 -I1 17 -m1 2.5 -M1 15.0 -A1 0.01 -d1 200 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 -tn 3 -t0 0.0 -t1 50 -t4 0 -tR 100 -xn 0 -xN 40.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 40 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H11 100.0 -H12 0.05 -Hp1 0.3 -Hp2 0.5 -HN 40 -tE 0.1 -G30 R -Fn 0

   PseudoData (output variables: -n 1 -v0 0): n_A at stationarity for a number of realizations (-tR 100); 

   MODEL: BEDDINGTON-DeANGELIS. Scanning 2D Parameter Space (Alpha_C_0, Nu_C_0):
   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 1 -v0 0 -G0 3 -G1 3 -sT 1.0E-06 -sN 300 -sP 2 -H9 2.5 -I0 16 -m0 0.5 -M0 5.0 -A0 0.01 -d0 200 -H10 10.0 -I1 17 -m1 2.5 -M1 15.0 -A1 0.01 -d1 200 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 -tn 3 -t0 0.0 -t1 50 -t4 0 -tR 100 -xn 0 -xN 40.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 40 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H11 100.0 -H12 0.05 -Hp1 0.3 -Hp2 0.5 -HN 40 -tE 0.1 -G30 R -Fn 0
   
   MODEL: BEDDINGTON-DeANGELIS. Scanning 2D Parameter Space (Chi_C_0, Eta_C_0):
   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 1 -v0 0 -G0 -3 -G1 3 -sT 1.0E-06 -sN 300 -sP 2 -H11 100.0 -I0 18 -m0 10 -M0 200.0 -A0 0.01 -d0 200 -H12 0.05 -I1 19 -m1 0.001 -M1 0.1 -A1 0.01 -d1 200 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 -tn 3 -t0 0.0 -t1 50 -t4 0 -tR 100 -xn 0 -xN 40.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 40 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H9 2.5 -H10 10.0  -Hp1 0.3 -Hp2 0.5 -HN 40 -tE 0.1 -G30 R -Fn 0
                                                      
   MODEL: BEDDINGTON-DeANGELIS. Scanning 2D Parameter Space (Alpha_C_0, Chi_C_0): 
   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 1 -v0 0 -G0 -3 -G1 3 -sT 1.0E-06 -sN 300 -sP 2 -H9 2.5 -I0 16 -m0 0.5 -M0 5.0 -A0 0.01 -d0 100 -H11 100.0 -I1 18 -m1 10 -M1 200.0 -A1 0.01 -d1 100 -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 -tn 3 -t0 0.0 -t1 50 -t4 0 -tR 100 -xn 0 -xN 40.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 40 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H10 10.0 -H12 0.05 -Hp1 0.3 -Hp2 0.5 -HN 40 -tE 0.1 -G30 R -Fn 0

   MODEL: HOLLING II. Scanning 2D Parameter Space (Alpha_C_0, Nu_C_0): 
   . ~$ ./DIFFUSION_HII_1D -y0 12 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 1 -v0 0 -G0 -3 -G1 3 -sT 1.0E-06 -sN 300 -sP 2 -H9 2.5 -I1 16 -m1 0.5 -M1 5.0 -A1 0.01 -d1 200 -H10 1.0 -I0 17 -m0 0.1 -M0 3.0 -A0 0.01 -d0 200  -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 -tn 3 -t0 0.0 -t1 1.5 -t4 0 -tR 100 -xn 0 -xN 20.0 -HN 20 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 15.0 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -Hp1 0.3750 -Hp2 0.5 -tE 0.1 -G30 R -Fn 0

   MODEL: HOLLING II nD. Scanning 2D Parameter Space (Alpha_C_0, Nu_C_0 from 1st Resource Type) 
   . ~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 5 -HM 1 -HX 1 -HY 1 \
                           -n 5 -v0 0 -v1 1 -v2 2 -v3 3 -v4 4 \
                           -G0 -3 -G1 3 \
                           -sT 1.0E-06 -sN 300 -sP 2 \
                           -H9 1.0 -I1 16 -m1 0.001 -M1 4.0 -A1 0.01 -d1 200 \
                           -H10 1.0 -I0 17 -m0 0.001 -M0 10.0 -A0 0.01 -d0 200 \
                           -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 \
                           -tn 2 -t0 0.0 -t1 0.8 -t4 0 -tR 100 -tE 0.1 -xn 0 -xN 20.0 -HN 20 \
                           -G2 1 -G3 0.0 -G4 4.0 -G5 1 -G6 0.0 -G7 10.0 \
                           -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -Hp1 0.3750 -Hp2 1.0 \
                           -G30 R -Fn 0

   . ~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 5 -HM 1 -HX 1 -HY 1 \
                           -n 5 -v0 0 -v1 1 -v2 2 -v3 3 -v4 4 \
                           -G0 2 -G1 2 \
                           -sT 1.0E-06 -sN 300 -sP 2 \
                           -H9 5.0 -I1 16 -m1 0.01 -M1 20.0 -A1 0.01 -d1 200 \
                           -H10 10.0 -I0 17 -m0 0.01 -M0 25.0 -A0 0.01 -d0 200 \
                           -iP 0 -en 0 -e0 426.012 -DP 0 -DC 0 -D0 0 -D1 0 -D2 0 -a0 0 \
                           -tn 2 -t0 0.0 -t1 0.25 -t4 0 -tR 10 -tE 0.1 -xn 0 -xN 20.0 -HN 20 \
                           -G2 1 -G3 0.0 -G4 15.0 -G5 1 -G6 0.0 -G7 10.0 \
                           -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -Hp1 1.0 -Hp2 1.0 \
                           -G30 R -Fn 0
   . Parameter Definitions: 
     -HuR -HuC are the jumping rates
     -H0 -H2 -H5  are the external immigration (Lambda_R_0, Lambda_R_1 and Lambda_C_0)
     -H20 is the establishment rate 
     -H1 -H3 -H6 are the death rates (Delta_R_0, Delta_R_1 for propagules, and Delta_C_0 for both searching and handling consumers)
     -H9 and -H10 are the Alpha_C_0 and Nu_C_0  Holling Type II model parameters 
     -H4 and -H17 are the production rates of propagules (Beta_R) and searching animals (Beta_C), respectively.  
     -xN is INITIAL_TOTAL_POPULATION (involved in fixing Ini Conditions: include.Initial_Conditions.[].c)
     -HN is No_of_INDIVIDUALS (involved in fixing some model parameters: see include.Parameter_Model.[].c )
     -xR  bool variable (0/1: Initial Conditions are not/yes re-scaled to the INITIAL_TOTAL_POPULATION value)
     -Fn 1/0 to generate pseudo data. Important: Use -Fn 0 to generate pseudo data. 
     -tn collects an input parameter controling the number time data points (usually two: Time_0 and Time_1). 
     -G30 R // Position of scale color bar: R, right side / L, left side / B, bottom side / T, top side .  
   
   . '#define' REPETITIONS (see below) tells the program how many times you will repeat the same experiment, 
      which consists in generating -tR [Realizations] from time -t0 [Initial Time] up to time -t1 [Last Time].
      However, the program may also ask you how many REPETITIONS to conduct at running time. 
*/
#define REPETITIONS 20000
// #define CONFIDENCE_INTERVALS

gsl_rng * r; /* Global generator defined in main.c */

void Writing_Empirical_Time_Vector(Parameter_Fitting * F, int n);
void Creating_HII_nD_Data_Matrix_from_Model ( Parameter_Table * , double **  );
void Creating_Standard_Data_Matrix_from_Model ( Parameter_Table * , double ** );
void Pointer_To_Function_Fitting_Structure (Parameter_Fitting * F, Parameter_Table * Table);
void Minimum_Parameter_2D_Scan(Parameter_Table * Table,
			       int No_of_POINTS_1, int Input_Parameter_1,
			       int No_of_POINTS_2, int Input_Parameter_2,
			       double * W_GRID,
			       double * Likelihood_Minimum,
			       double * x_Value, 
			       double * y_Value); 
void Profiling_2D_Scanned_Function(Parameter_Table * Table,
				   int No_of_POINTS_1, int Input_Parameter_1,
				   int No_of_POINTS_2, int Input_Parameter_2,
				   double * W_GRID,
				   double Initial_Parameter_Val_1, 
				   double Initial_Parameter_Val_2, 
				   double * Profile_x_Data, double * x_Data, 
				   double * Profile_y_Data, double * y_Data);
void Profiling_2D_Scanned_Function_Maximum_Minimum(Parameter_Table * Table,
						   int No_of_POINTS_1, int Input_Parameter_1,
						   int No_of_POINTS_2, int Input_Parameter_2,
						   double * W_GRID,
						   double * Profile_x_Data_B, double * x_Data, 
						   double * Profile_y_Data_L, double * y_Data,
						   double * Profile_x_Data_T, 
						   double * Profile_y_Data_R);
void Confidence_Interval_Parameter_Ratio( Parameter_Table * Table,
					  int Input_Parameter_0,
					  int Input_Parameter_1, 
					  double x_B, double * x_CI_B,
					  double x_T, double * x_CI_T,
					  double y_L, double * y_CI_L,
					  double y_R, double * y_CI_R,
					  double * Parameter_Ratio, 
					  double * Parameter_Ratio_CI ); 

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

float * customized_contour_levels_1( Parameter_CPGPLOT * C, double Min_Value)
{
    int i;

    /* Two contour levels */
    C->NC = 2;
    float * clevels = (float *)calloc( C->NC, sizeof(float) );
    /* clevels[0] = Min_Value + 2.0;   */
    /* clevels[1] = Min_Value + 100.0; */
    clevels[0] = Min_Value + 0.001;
    clevels[1] = Min_Value + 2.0;
    
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
  double Likelihood_Minimum, x_Val, y_Val;

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
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C( &Table );
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( &Table, Index_Output_Variables );
  printf(" Parameter_Table structure has been correctly allocated and initiated\n");

  Parameter_Model * City_Par_Values = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (&Table, City_Par_Values);
  printf(" Parameter_Model structure 'City_Par_Values' has been correctly allocated and initiated\n");
  
  /* B E G I N : Reserving memmory for Parameter Space */
  Parameter_Space * Error_Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  assert( No_of_ERROR_PARAMETERS == 0 );
  Table.E_Space = Error_Space;
 
  Parameter_Space * Initial_Condition_Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
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
    Empirical_Data_Matrix[i] = (double *)calloc( Realizations, sizeof(double) );
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
    					                                TYPE_of_TIME_DEPENDENCE,                // -t4 1
    					                                TYPE_0_PARAMETERS,                      // -D0 0
    					                                TYPE_1_PARAMETERS,                      // -D1 1
    					                                TYPE_2_PARAMETERS,                      // -D2 0
    					                                No_of_COVARIATES,                       // -DC 0
    					                                dependent_parameter, forcing_pattern,   // -P0 16
    					                                Matrix_of_COVARIATES,
    					                                Type_1_Parameter_Values,
    					                                Time_Empirical_Vector);
  }
  /*     E N D -----------------------------------------------------------------------*/

  if(Table.TYPE_of_MODEL == 12 || Table.TYPE_of_MODEL == 13 || Table.TYPE_of_MODEL == 14) 
    /* Models where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);

  if(Table.TYPE_of_MODEL == 16) { // DIFFUSION_HIIl_nD
    /* A model where the TOTAL_No_of_CONSUMERS is a CONSTANT and consumers       */
    /* feed on multiple resources: extra parameters per resource type are needed */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
    Resetting_Alpha_Nu_Vectors_Constant (&Table); /* Nu's and Alpha's are all equal */
    Resetting_Multiresource_Levels (&Table);      /* Creating Vector of Thetas      */
    Writing_Alpha_Nu_Theta_Vectors(&Table);  
  }
  
  /* B E G I N : Observed Data Control Initization (Reserving Memmory)          */
  Observed_Data * Data = (Observed_Data *)calloc(1, sizeof(Observed_Data));
  Observed_Data_Alloc( Data, SUB_OUTPUT_VARIABLES, Time.Realizations);
  printf(" Observed_Data structure has been correctly allocated\n");
  /*     E N D : -------------------------------------------------------------- */
  /* B E G I N : Reserving memmory for Parameter Fitting Structure ------------ */
  Parameter_Fitting * F = (Parameter_Fitting*)calloc(1,sizeof(Parameter_Fitting));
  F->Data                  = Data;
  F->Space                 = Space;
  F->Table                 = &Table;
  F->Minimization          = 0;     
  F->Bounded_Parameter_Set = 1;
  Pointer_To_Function_Fitting_Structure (F, &Table);
  // For instance (see function definition below),  
  // F->Function              = GSL_Function_to_Minimize_Binomial_Free_Consumers; 
  /* double Function_to_Minimize( Parameter_Table * Table ) 
     This will be the function that will be past to the scanning function
     void generic_Function_Parameter_2Dim_Scan_Improved(Parameter_Table * Table, ...)
     The function 'Function_to_Minimize' is just a correctly prototyped wrapper to 
     generalically scan the function that the pointer 'F->Function' points to 
     (see Funtion_to_Miminize definintion in Optimization Library in 
     ./Library directory). 
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
  /*     E N D : ----------------------------------------------------------------- */
  
  int No_of_COLS;
  if( No_of_FILES > 0) No_of_COLS = F_y_GRID[0];        /* -Y0 [VALUE] */
  else                 No_of_COLS = Realizations;       /* -tn [VALUE] */

  Initial_Value_0 = AssignStructValue_to_VectorEntry(Input_Parameter_1, &Table);
  Initial_Value_1 = AssignStructValue_to_VectorEntry(Input_Parameter_2, &Table);  

  FILE * FP_x = fopen("Confidence_Intervals_0_HII.dat", "w"); 
  FILE * FP_y = fopen("Confidence_Intervals_1_HII.dat", "w"); 

  FILE * FP_mle_x = fopen("MLE_0_HII.dat", "w"); 
  FILE * FP_mle_y = fopen("MLE_1_HII.dat", "w");

  double * x_Val_mle = (double *)calloc(REPETITIONS, sizeof(double)); 
  double * y_Val_mle = (double *)calloc(REPETITIONS, sizeof(double));

  int No_of_REPETITIONS = REPETITIONS;  

  printf("Enter and integer-valued number, REPETITIONS, please... ");
  printf("It will be the number of repeated MLE estimation procedures conducted.\n");
  scanf("%d", &No_of_REPETITIONS); 

  for(k=0; k<No_of_REPETITIONS; k++) {
    /* Every new repetition requires setting up Time_Vector_Real back to its initial values
       because, after each repetition, this vector has taken the values of real final times
       (Time_Current values of the last step) for every realization. 
    */

    for(i = 0; i<Time.Realizations; i++) {
      Time.Time_Vector_Real[i][0] = Time_0;
      Time.Time_Vector_Real[i][1] = Time_0 + (double)(i+1) * (Time_1 - Time_0)/(double)(Time.Realizations);
    }
    
    /* 
      Time_0 = 0.0; 
      Time_1 = Table.Tiempo_Propio;  
      for(i = 0; i<Time.Realizations/2; i++) {
        Time.Time_Vector_Real[i][0] = Time_0;
        Time.Time_Vector_Real[i][1] = Time_0 + (double)(i+1) * (Time_1 - Time_0)/(double)(Time.Realizations);
      }

      Time_0 = 9.0 * Table.Tiempo_Propio; 
      Time_1 = 10.0 * Table.Tiempo_Propio; 
      for(i = Time.Realizations/2; i<Time.Realizations; i++) {
        Time.Time_Vector_Real[i][0] = Time_0;
        Time.Time_Vector_Real[i][1] = Time_0 + (double)(i+1) * (Time_1 - Time_0)/(double)(Time.Realizations);
      }
    */

    printf("\t Total number of experiment repetitions: %d (current repetion: %d)\n",
	          No_of_REPETITIONS, k);
    printf("\t %d stochastic realizations are sampled from t_0=%.2g to t_1\n",
	          Realizations, Time_0);
    
    for(i = 0; i<Time.Realizations; i++) 
      printf("\n t_1 (Relization %d) = %1.3g\n", i, Time.Time_Vector_Real[i][1]);
    
    printf("\n For each realization, a different final time is used!!!\n"); 
    
    /* In order to produce pseudata always with the same true parameters, they need to be set up again. 
       because the 2D scan changes its values in the Parameter_Table structure  
    */
    AssignVectorEntry_to_Structure(&Table, Input_Parameter_1, Initial_Value_0);
    AssignVectorEntry_to_Structure(&Table, Input_Parameter_2, Initial_Value_1);

    /* As well as anything that depends on them (such as the Nu and Alpha vectors )   */
    if(Table.TYPE_of_MODEL == 16) { // DIFFUSION_HIIl_nD
      /* A model where the TOTAL_No_of_CONSUMERS is a CONSTANT and consumers          */
      /* feed on multiple resources: extra parameters per resource type are needed    */
      Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
      Resetting_Alpha_Nu_Vectors_Constant (&Table); /* Nu's and Alpha's are all equal */
      Resetting_Multiresource_Levels (&Table);      /* Creating Vector of Thetas      */
      Writing_Alpha_Nu_Theta_Vectors(&Table);  
    }
    /* B E G I N : Generation of the matrix of pseudo-data (always with the same true values)  */
    if( No_of_FILES == 0) M_O_D_E_L___S_T_O( &Table ); 
    /* E N D     : --------------------------------------------------------------------------- */

    if ( No_of_FILES > 0 )
      Reading_Standard_Data_Matrix_from_File( OBSERVED_DATA_FILE,
					                                    Empirical_Data_Matrix,
					                                    &SUB_OUTPUT_VARIABLES,
					                                    No_of_COLS, 
					                                    0, Name_of_Rows,
					                                    1, Time.Time_Vector );
    else
      /* Storing Pseudo-data in Empirical Data Matrix   */
      /* (this process is model specific, in general)   */
      Creating_HII_nD_Data_Matrix_from_Model ( &Table, Empirical_Data_Matrix );
      // Creating_Standard_Data_Matrix_from_Model ( &Table, Empirical_Data_Matrix );

    Realizations = Table.T->Realizations; /* In MODEL DIFFUSION_HII_nD, we have 
                                             as many realizations as different 
                                             final times 
                                          */
    double * Realizations_Vector = (double *)calloc(Realizations, sizeof(double));
    for (i=0; i<Realizations; i++)
      Realizations_Vector[i] = (double)i + 1.0;
    
    // #if defined VERBOSE
    Writing_Empirical_Time_Vector(F, Realizations);
    Writing_Standard_Data_Matrix( Empirical_Data_Matrix,
				                          SUB_OUTPUT_VARIABLES, Realizations,
				                          1, Name_of_Rows,
				                          0, Realizations_Vector);
    // #endif

    /* B E G I N : Observed Data Control Initization (Initializing value)       */
    Observed_Data_Initialization( Data, SUB_OUTPUT_VARIABLES,
				                          Realizations, Empirical_Data_Matrix,
				                          "" );
    printf(" Observed_Data structure has been correctly initiated with an instance\n");
    printf(" of (real or pseudo) empirical data \n");
    Print_Press_Key(0, 0, ".");
    /*     E N D : -------------------------------------------------------------- */

    /* Back to input argument values in Table */
    AssignVectorEntry_to_Structure(&Table, Input_Parameter_1, Initial_Value_0);
    AssignVectorEntry_to_Structure(&Table, Input_Parameter_2, Initial_Value_1);
    P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( &Table, City_Par_Values) ;

    /* B E G I N : Main Function Call (the core of this code) -------------------------*/
    double * W_GRID = (double *)malloc( No_of_POINTS_1 * No_of_POINTS_2 * sizeof(double) );
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
    Minimum_Parameter_2D_Scan(&Table,
			                        No_of_POINTS_1, Input_Parameter_1,
			                        No_of_POINTS_2, Input_Parameter_2,
			                        W_GRID,
			                        &Likelihood_Minimum, &x_Val, &y_Val);

    fprintf(FP_mle_x, "%g\t%g\n", Initial_Value_0, x_Val); 
    fprintf(FP_mle_y, "%g\t%g\n", Initial_Value_1, y_Val); 

    x_Val_mle[k] = x_Val; 
    y_Val_mle[k] = y_Val;

    printf(" Output values of Minimum Parameter 2D Scan(...) function:\n");
    printf(" Optimal Negative logLikelihood: %g\n", Likelihood_Minimum); 
    printf("   %s=%f\t", Table.Symbol_Parameters[Input_Parameter_1], x_Val);
    printf("   %s=%f\n", Table.Symbol_Parameters[Input_Parameter_2], y_Val);
    printf(" (This is the mininum loglikelihood value over the 2D scan\n");
    printf("  within the subparameter space analyzed)\n");
    /*   E N D : Main Function Call -------------------------------------------------*/
    
  if(No_of_REPETITIONS < 10) {
   /* B E G I N :   Drawing heat maps ------------------------------------------------*/    
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
    /* B E G I N : Annotating the countours by hand and drawing true value */
    // cpgptxt(float x, float y, float angle, float fjust,  const char *text);

    cpgslw(2); 
    /* cpgptxt(0.00002, 62.0, 0.0, 0.0,  "1.0"); */
    /* cpgptxt(0.0001, 78.0, 0.0, 0.0,   "2.5"); */
    /* cpgptxt(0.0003, 92.0, 0.0, 0.0,  "5.0");  */
    
    float x_Value, y_Value; 
    float * xs; 
    float * ys; 

    x_Value = Parameter_Model_into_Vector_Entry(Input_Parameter_1, City_Par_Values);
    y_Value = Parameter_Model_into_Vector_Entry(Input_Parameter_2, City_Par_Values);
    
    cpgslw(3);  /* Line width changing to 3     */
    cpgsci(2);  /* Color Index changing to 12   */
    cpgpt1(x_Value, y_Value, 23);  /* Symbol 23 */
  
    xs = (float *)calloc(2, sizeof(float) );
    ys = (float *)calloc(2, sizeof(float) );
    xs[0] = 0.5* x_Value;  xs[1] = x_Value; /* A 40 % reduction */ 
    ys[0] = y_Value;       ys[1] = y_Value;
    
    cpg_XY_same_arrow( 2, xs, ys, 5, 1, 4);

    printf("%s = %g[ %f ]\t", Table.Symbol_Parameters[Input_Parameter_1], x_Val, x_Value);
    printf("%s = %g[ %f ]\n", Table.Symbol_Parameters[Input_Parameter_2], y_Val, y_Value);
    
    x_Value = (float)x_Val;
    y_Value = (float)y_Val; 
    
    cpgslw(3);  /* Line width changing to 3     */
    cpgsci(2);  /* Color Index changing to 12   */
    cpgpt1(x_Value, y_Value, 23);  /* Symbol 23 */
  
    xs = (float *)calloc(2, sizeof(float) );
    ys = (float *)calloc(2, sizeof(float) );
    xs[0] = 1.40 * x_Value;  xs[1] = x_Value; /* A 40 % reduction */ 
    ys[0] = y_Value;         ys[1] = y_Value;
    
    cpg_XY_same_arrow( 2, xs, ys, 5, 1, 4);
    // cpg_XY_same_arrow( N, xs, ys, CPG->color_Index, CPG->type_of_Line, CPG->type_of_Width );
    
    free(xs);
    free(ys);
    Print_Press_Key(1,0,".");
    /* E N D : Annotating the countours by hand and drawing true value */
  /* E N D :  Drawing heat maps ----------------------------------------------------------*/
  /* B E G I N :   Drawing likelihood profiles -------------------------------------------*/  
  #if defined CONFIDENCE_INTERVALS
    double * x_Data  = (double *)calloc(No_of_POINTS_1, sizeof(double) ); 
    double * y_Data  = (double *)calloc(No_of_POINTS_2, sizeof(double) );
    double * Profile_x_Data  = (double *)calloc(No_of_POINTS_1, sizeof(double) ); 
    double * Profile_y_Data  = (double *)calloc(No_of_POINTS_2, sizeof(double) );
   
    Profiling_2D_Scanned_Function(&Table,
    				                      No_of_POINTS_1, Input_Parameter_1,
    				                      No_of_POINTS_2, Input_Parameter_2,
    				                      W_GRID,
    				                      Initial_Value_0,
    				                      Initial_Value_1,
    				                      Profile_x_Data, x_Data,
    				                      Profile_y_Data, y_Data);
    int SAME_PLOT = 0;
    int SCALE_X   = 0;
    int SCALE_Y   = 1;
    Table.CPG->color_Index   = 2;
    Table.CPG->type_of_Line  = 1;
    Table.CPG->type_of_Width = 5;
    Table.CPG->type_of_Symbol= 1;
    Table.CPG->CPG_RANGE_Y_0 = Likelihood_Minimum - 10.0;
    Table.CPG->CPG_RANGE_Y_1 = Likelihood_Minimum + 15.0;
    
    CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T(Table.CPG,
    							                                      SAME_PLOT,
    							                                      No_of_POINTS_1,
    							                                      x_Data, Profile_x_Data,
    							                                      Table.CPG->X_label,
    							                                      "Negative Log Likelihood",
    							                                      "",
    							                                      SCALE_X, SCALE_Y );
    Print_Press_Key(1,0,".");
    CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T(Table.CPG,
    							                                      SAME_PLOT,
    							                                      No_of_POINTS_2,
    							                                      y_Data, Profile_y_Data,
    							                                      Table.CPG->Y_label,
    							                                      "Negative Log Likelihood",
    							                                      "",
    							                                      SCALE_X, SCALE_Y );
    Print_Press_Key(1,0,".");
    /*    E N D :   Drawing likelihood profiles -------------------------------------------*/
    
    printf("Calculating Confidence Intervals (Input_Parameter_1 and Input_Parameter_2)...\n");
    Print_Press_Key(1,0,".");
    /* B E G I N :  Calculating Confidence Intervals (Input_Parameter_1 and Input_Parameter_2) */
    printf("Confidence Intervals from Likelihood Profile:\n");
    double x_MLE = 0.0;
    double LIKELIHOOD_JUMP = 2.0;
    double * x_CI = (double *)calloc(2, sizeof(double) ); 
    Confidence_Intervals_from_Likelihood_Profile( Profile_x_Data, x_Data, No_of_POINTS_1,
						                                      LIKELIHOOD_JUMP, 
						                                      &x_MLE, x_CI );
    double y_MLE = 0.0;
    double * y_CI = (double *)calloc(2, sizeof(double) ); 
    Confidence_Intervals_from_Likelihood_Profile( Profile_y_Data, y_Data, No_of_POINTS_2,
						                                      LIKELIHOOD_JUMP, 
						                                      &y_MLE, y_CI );
    /* True Value,  MLE,  [ CI[0],  CI[1] ] */
    fprintf(FP_x, "%g\t%g\t[%g, %g]\n", Initial_Value_0, x_MLE, x_CI[0], x_CI[1]); 
    fprintf(FP_y, "%g\t%g\t[%g, %g]\n", Initial_Value_1, y_MLE, y_CI[0], y_CI[1]); 
    
    printf("[True Value: %s = %g]\t%s=%f\t[%f, %f]\n",
	  Table.Symbol_Parameters[Input_Parameter_1], Initial_Value_0,
	  Table.Symbol_Parameters[Input_Parameter_1], x_MLE, x_CI[0], x_CI[1]);
    
    printf("[True Value: %s = %g]\t%s=%f\t[%f, %f]\n",
	  Table.Symbol_Parameters[Input_Parameter_2], Initial_Value_1,
	  Table.Symbol_Parameters[Input_Parameter_2], y_MLE, y_CI[0], y_CI[1]);
    
    free(x_CI); free(y_CI);    
    free(x_Data);         free(y_Data);
    free(Profile_x_Data); free(Profile_y_Data);
    /*     E N D :  Calculating Confidence Intervals ------------------------------------------*/
  #endif
    Print_Press_Key (1, 1, "Next Experiment Repetition");
  }
#endif
    
  free (W_GRID);
  free (Realizations_Vector); 
}

if(No_of_REPETITIONS > 500) {  
  Table.CPG->CPG_RANGE_X_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_min );
  Table.CPG->CPG_RANGE_X_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_MAX );
  int No_of_BINS = No_of_REPETITIONS/100;
  int SAME = 0; int BAR = 1; 
  CPGPLOT___X_Y___H_I_S_T_O_G_R_A_M___N_O_R_M_A_L_I_Z_E_D( SAME,
                                                           Table.CPG, 
                                                           No_of_REPETITIONS, x_Val_mle, No_of_BINS, 
                                                           Table.Symbol_Parameters[Input_Parameter_1], 
                                                           "Frequency", "", 
                                                           0, 0, 
                                                           BAR ); 
  Print_Press_Key(1, 1, "Next Normalized Histrogram"); 

  Table.CPG->CPG_RANGE_X_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_min );
  Table.CPG->CPG_RANGE_X_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_MAX ); 
  CPGPLOT___X_Y___H_I_S_T_O_G_R_A_M___N_O_R_M_A_L_I_Z_E_D( SAME,
                                                           Table.CPG, 
                                                           No_of_REPETITIONS, y_Val_mle, No_of_BINS, 
                                                           Table.Symbol_Parameters[Input_Parameter_2], 
                                                           "Frequency", "", 
                                                           0, 0,
                                                           BAR ); 
}

  Print_Press_Key(1, 1, "Freeing All Memmory...");
  free(x_Val_mle);  free(y_Val_mle);
  fclose(FP_x);     fclose(FP_y);
  fclose(FP_mle_x); fclose(FP_mle_y);

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
  for (i=0; i<SUB_OUTPUT_VARIABLES; i++)  
    free(Empirical_Data_Matrix[i]);
  free(Empirical_Data_Matrix);

  free(TIME_PARAMETERS_FILE);
  free(OBSERVED_DATA_FILE);
  free(PARAMETER_SET_FILE);

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */
  
  printf("\nEnd of progam\n");
  return (0);
}

void Writing_Empirical_Time_Vector(Parameter_Fitting * F, int n)
{
  int i;

  for(i=0; i<n; i++)
    printf("T[%d]  ", i);
  printf("\n");

  for(i=0; i<n; i++)
    printf("%g  ", F->Data->Time_Data_Vector[i]);  
  printf("\n");
} 

void Minimum_Parameter_2D_Scan (Parameter_Table * Table,
			                          int No_of_POINTS_1, int Input_Parameter_1,
			                          int No_of_POINTS_2, int Input_Parameter_2,
			                          double * W_GRID,
			                          double * Likelihood_Minimum,
			                          double * x_Value_MLE, 
			                          double * y_Value_MLE)
{
  /* This function looks for a minimum of the function over the 2D subparameter space that 
     has been scanned (as defined in Table->S Parameter_Space structure
  */ 
  int n, k, j, i;
  int k_MIN, j_MIN; 
  double Minimum_Value, Value, Value_0, Value_1;
	
  Parameter_Space * S = Table->S;
	/* BEGIN : Allocating memory for a 2D array to save negative loglikelioods */  
	double      ** z_SOL  = (double **)malloc( No_of_POINTS_2 * sizeof(double *) );
	for( i = 0; i < No_of_POINTS_2; i++){
	  z_SOL[i] = (double *)malloc( No_of_POINTS_1 * sizeof(double) );
	}
	double * x_Data  = (double *)malloc(No_of_POINTS_1 * sizeof(double) ); 
	double * y_Data  = (double *)malloc(No_of_POINTS_2 * sizeof(double) ); 
	/*   END : --------------------------------------------------------------- */

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

  #if defined VERBOSE
	printf("Optimal Negative logLikelihood: %g\n", Minimum_Value); 
	printf("%s=%f\t", Table->Symbol_Parameters[Input_Parameter_1], x_Data[j_MIN]);
	printf("%s=%f\n", Table->Symbol_Parameters[Input_Parameter_2], y_Data[k_MIN]);
	#endif

	* Likelihood_Minimum = Minimum_Value; 
	* x_Value_MLE = x_Data[j_MIN];
	* y_Value_MLE = y_Data[k_MIN];
	
	/* BEGIN : Free memory!!!                                     * * * * * * */
	for( i = 0; i < No_of_POINTS_2; i++) 
    free(z_SOL[i]); 
	free(z_SOL);
	
  free(x_Data); free(y_Data); 
	/*   END : Free memmory!!!                                    * * * * * * */
}

void Profiling_2D_Scanned_Function( Parameter_Table * Table,
				                            int No_of_POINTS_1, int Input_Parameter_1,
				                            int No_of_POINTS_2, int Input_Parameter_2,
				                            double * W_GRID,
				                            double Initial_Parameter_Val_1, 
				                            double Initial_Parameter_Val_2, 
				                            double * Profile_x_Data, double * x_Data, 
				                            double * Profile_y_Data, double * y_Data)
{
	int n, k, j, i;
	int k_VAL, j_VAL, k_BOOL, j_BOOL; 
	double Value, Value_0, Value_1, Step_y, Step_x;
	
	Parameter_Space * S = Table->S;
	/* BEGIN : Allocating memory for saving data to plot a bifurcation  * * * * * * */
	/*         diagram for each variable  * * * * * * * * * * * * * * * * * * * * * */  
	double      ** z_SOL  = (double **)malloc( No_of_POINTS_2 * sizeof(double *) );
	for( i = 0; i < No_of_POINTS_2; i++)
	  z_SOL[i] = (double *)malloc( No_of_POINTS_1 * sizeof(double) );
	/*   END : Allocating memory for saving dynamical data * * * * * * * * * * * *  */

	k_BOOL = 0; j_BOOL = 0; 
	n = 0; 
	for( k = 0; k < No_of_POINTS_2; k++ ) {

	  Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min );
	  Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX );

	  Step_y  = (Value_1 - Value_0)/(double)(No_of_POINTS_2 - 1);
     
	  Value = Value_0 + k * (Value_1 - Value_0)/(double)(No_of_POINTS_2 - 1);
	  y_Data[k]= Value;

	  if( Initial_Parameter_Val_2 >= Value && Initial_Parameter_Val_2 < Value + Step_y ){
	    k_VAL = k;
	    k_BOOL = 1; 
	  }
	    
	  Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min );
	  Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX );

	  Step_x  = (Value_1 - Value_0)/(double)(No_of_POINTS_1 - 1);
	  
	  for( j = 0; j < No_of_POINTS_1; j++ ){
	    
	    Value = Value_0 + j * (Value_1 - Value_0)/(double)(No_of_POINTS_1 - 1);
	    
	    x_Data[j] = Value;

	    z_SOL[k][j]    = W_GRID[n++]; 
	    
	    if(Initial_Parameter_Val_1 >= Value && Initial_Parameter_Val_1 < Value + Step_x){
	      j_VAL = j;
	      j_BOOL = 1; 
	    }
	  }
	}

	if (j_BOOL == 1 && k_BOOL == 1 ) {
	  for( k = 0; k < No_of_POINTS_2; k++ ) Profile_y_Data[k] = z_SOL[k][j_VAL]; 
	  
	  for( j = 0; j < No_of_POINTS_1; j++ ) Profile_x_Data[j] = z_SOL[k_VAL][j];
	  
	  printf("Parameter values from the discretized profiles:\n"); 
	  printf("%s=%f\t", Table->Symbol_Parameters[Input_Parameter_1], x_Data[j_VAL]);
	  printf("%s=%f\n", Table->Symbol_Parameters[Input_Parameter_2], y_Data[k_VAL]);	
	  printf("Initial Parameter values (as the input arguments to this function):\n"); 
	  printf("%s=%f\t", Table->Symbol_Parameters[Input_Parameter_1], Initial_Parameter_Val_1);
	  printf("%s=%f\n", Table->Symbol_Parameters[Input_Parameter_2], Initial_Parameter_Val_2);	
	}
	else {
	  printf(" Initial parameter values are out of range\n");
	  printf(" The program will exit\n");
	  exit(0); 
	}
	
	/* BEGIN : Free memory!!!                                     * * * * * * */
	for( i = 0; i < No_of_POINTS_2; i++) free(z_SOL[i]); 
	free(z_SOL);
	/*   END : Free memmory!!!                                    * * * * * */
}
      
void Creating_Standard_Data_Matrix_from_Model ( Parameter_Table * Table,
						                                    double **  Empirical_Data_Matrix )
{
  /* This is just to save the output variables corresponding to the i-th realization 
     generated by the function  M_O_D_EL___S_T_O( &Table ) into a straightforward 
     Empirical Data Matrix containing the 'Psedo Data' as have been generated in the 
     i-th realization. 
  */
  int k,i;
  int j    = Table->T->I_Time - 1; /* Value corresponding to very last time */
  int F;       /* F: Number of errors at the j-th time across realizations; */
  int count, i_valid, valid_realizations;
  int No_of_POINTS; 
  
  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++) {

    count   = 0;
    i_valid = 0;
    F       =  Table->T->Realizations - Table->T->count[j];

    for(i=0; i<Table->T->Realizations; i++) {
      
      if( Table->T->Variable[i][k][j] == 0.0 ) {
	      if (count < F)
	        count++;
	      else
	        Empirical_Data_Matrix[k][i_valid++] = Table->T->Variable[i][k][j];
      }
      else 
	      Empirical_Data_Matrix[k][i_valid++] = Table->T->Variable[i][k][j];
    }
    valid_realizations = i_valid;
  }

  Table->T->Realizations = valid_realizations; 
}

void Creating_HII_nD_Data_Matrix_from_Model ( Parameter_Table * Table,
						                                  double **  Empirical_Data_Matrix )
{
  /* 
     This function saves the output variables corresponding to the realizations 
     generated by the function  M_O_D_EL___S_T_O( &Table ) into a straightforward 
     Empirical Data Matrix containing just the 'Pseudo Data' required for the 
     particular likelihood of the MODEL DIFFUSION_HII_nD to be calculated. 
     This likelihood procedure uses only two times, I_Time = 2 (-tn 2 in argument list)
  */
  double t, p;
  int k, i;
  int No_of_POINTS; 

  Parameter_Fitting * F = (Parameter_Fitting *)Table->Fitting_Data; 

  assert(Table->T->I_Time == 2);

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++) {

    for(i=0; i<Table->T->Realizations; i++) {
      /* For each time, the i-th realization corresponding to the i-th time 
         is stored in Empirical Data Matrix 
      */ 
	    Empirical_Data_Matrix[k][i] = Table->T->Variable[i][k][1];

      /* Associated times */
      F->Data->Time_Data_Vector[i] = Table->T->Time_Vector_Real[i][1];          
    }
  }

  /* Estimating Average Values at time t: <n> = p_i(t) * No_of_CONSUMERS */
  double ** n_AVE = (double **)calloc(Table->No_of_RESOURCES, sizeof(double *));
  for(k = 0; k<Table->No_of_RESOURCES; k++)
    n_AVE[k] = (double *)calloc(Table->T->Realizations, sizeof(double)); 
  
  double Theta_Sum, Sum_ThetaNu_Ratio;
  Theta_Nu_Sums (Table->Theta_Consumers, Table->Nu_Consumers, 
                 &Theta_Sum, &Sum_ThetaNu_Ratio, 
                 Table->No_of_RESOURCES);

  for(k = 0; k<Table->No_of_RESOURCES; k++)
    for(i = 0; i<Table->T->Realizations; i++) {
      t = F->Data->Time_Data_Vector[i];

      p = Function___poft (Table->Theta_Consumers[k], Table->Nu_Consumers[k], 
                           Sum_ThetaNu_Ratio, Theta_Sum, 
                           t);
      n_AVE[k][i] = p*Table->TOTAL_No_of_CONSUMERS;
    }

  /* Saving pseudo-data matrix in file Empirical_Data_Matrix.dat */
  FILE * fp = fopen("Empirical_Data_Matrix.dat", "w");
  for(i=0; i<Table->T->Realizations; i++) 
    fprintf(fp, "%g\t", F->Data->Time_Data_Vector[i]);
  
  fprintf(fp, "\n");

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++) {
      for(i=0; i<Table->T->Realizations; i++) {
        fprintf(fp, "%g\t", Empirical_Data_Matrix[k][i]);
      }
    fprintf(fp, "\n");
  }
  fclose(fp);

  /* Saving theoretical matrix (with expected values) in file Theoretical_Expected_Matrix.dat */
  fp = fopen("Theoretical_Expected_Matrix.dat", "w");
  for(i=0; i<Table->T->Realizations; i++) 
    fprintf(fp, "%g\t", F->Data->Time_Data_Vector[i]);

  fprintf(fp, "\n");

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++) {
      for(i=0; i<Table->T->Realizations; i++) {
        fprintf(fp, "%g\t", n_AVE[k][i]);
      }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void Pointer_To_Function_Fitting_Structure (Parameter_Fitting * F, Parameter_Table * Table)
{
  /* This function assigns the correct pointer to the Fitting Structure
     function member in agreement with the model at work. 
     Only a couple of models have these functions defined. These are 
     analytical expression for the distributions at stationarity and/or
     at the system evolves in time.
  */
  switch( Table->TYPE_of_MODEL ) {

  case 9: /* DIFFUSION_HII_2D * * * * * * * * * * * * * * * * * * * * * * */   
    F->Function = GSL_Function_to_Minimize_Binomial_Free_Consumers; 
    break;
    
  case 12: /* DIFFUSION_HII_1D * * * * * * * * * * * * * * * * * * * * * * */
    F->Function = GSL_Function_to_Minimize_Binomial_Free_Consumers;
    break;

  case 13: /* DIFFUSION_BD_2D * * * * * * * * * * * * * * * * * * * * * * */
    if(Table->SUB_OUTPUT_VARIABLES == 1)
      F->Function = GSL_Function_to_Minimize_Beddington_DeAngelis_Marginal_0;
    else {
      assert(Table->SUB_OUTPUT_VARIABLES == 2); 
      F->Function = GSL_Function_to_Minimize_Beddington_DeAngelis;
    }
    break;
  
  case 16: /* DIFFUSION_HII_nD * * * * * * * * * * * * * * * * * * * * * * */
    F->Function = GSL_Function_to_Minimize_Multinomial_Free_Consumers;
    break;
      
  default:
    printf(" This TYPE_of_MODEL (%d) code does not have\n", Table->TYPE_of_MODEL);
    printf(" a defined likelihood to fit data\n"); 
    printf(" Ony likelihoods for models (9, 12, 13, 16) are correctly defined\n");
    printf(" (DIFFUSION HII 2D, DIFFUSION HII 1D, DIFFUSION BD 2D, and DIFFUSION_HII_nD)\n");
    printf(" Check input argument list!!!\n");
    exit(0);
  } 
}  
