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

   If the parametric function is a likelihood, it also depends on some data and 
   a Observed_Data structure should be defined and a corresponding data file read. 
   See an example of this in: 

   ~/PROJECT_STOCHASTIC_DIFFUSION/MODEL_CALCULATIONS/VISUAL_2DIM_PARAMETER_SCAN/main.c

   Otherwise, if the parametric function depends only on the parametric configuration, 
   the main code becomes simpler, as done here. Examples of such funcitons are 
   R_0 or any output variable that depends on the stationary state of the system. 

   In this example, for the DIFFUSION_1R1C model, a three-valued parametric function 
   is defined representing the qualitative different dynamics obtained for different
   parameteric configurations (stable end points, damped oscillatins, extinction of 
   consumers...): 

   . Function (Parameter_Table * Table) == 0, for only resources

   . Function (Parameter_Table * Table) == 1, for stable coexistce of resources
   and consumers

   . Function (Parameter_Table * Table) == 2, for stable coexistce of resources
   and consumers and damped oscillations 

   Other functions, for other models, can be defined in a similar way. 

   Compilation:
   
   . ~$ make MODEL=[TYPE_of_MODEL] 
   
   Execution Examples:
Index  Argument   Parameter Definition
-----  --------   --------------------
 0:   -HuR 0.0    Local Diffusion of Resources                
20:   -HuC 0.0    Local Diffusion of Consumers  (also -H13)
 5:   -HS  1      Number of Resources 
 6:   -H0 0.0     External Immigration of Resources
12:   -H5 0.0     External Immigration of Consumers
10:   -HK 2000.0  Resource Carrying Capacity
11:   -H4 1.5     Resource Growth Rate
 7:   -H1 0.5     Resource Death Rate 
13:   -H6 1.0     Consumer Death Rate
16:   -H9 5.0     Consumer Attack Rate
17:   -H10 1.0    Consumer Handling Time
18:   -H11 0.0    Rate of Triplet Formation 
19:   -H12 1.0    Rate of Triplet Degradation

   Extension of MacArthur and Rosenweig (4D: R, A, RA, ARA) with consumer interference. Notice -H11 [Xhi] and -H12 [Eta]. If these are zero, the model collapses into a 3D without consumer interference
   . ~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -I0 16 -m0 2.0 -M0 15.0 -A0 0.01 -d0 100  -I1 17 -m1 1.0 -M1 5.0 -A1 0.01 -d1 100 -iP 0 -en 0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 1.5 -H1 0.5 -H6 5.0 -H9 10.0 -H10 2.0 -H11 0.0 -H12 0.0
   
   MacArthur and Rosenzweig (3D: R, A, RA). Holling Type II
   . ~$ ./DIFFUSION_MR -y0 7 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -I0 16 -m0 2.0 -M0 15.0 -A0 0.01 -d0 100  -I1 17 -m1 1.0 -M1 5.0 -A1 0.01 -d1 100 -iP 0 -en 0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 1.5 -H1 0.5 -H6 5.0 -H9 10.0 -H10 2.0 

   MacArthur and Rosenzweig (3D: R, A, RA). Holling Type II
   . ~$ ./DIFFUSION_MR -y0 7 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -I0 16 -m0 0.5 -M0 20.0 -A0 0.01 -d0 200  -I1 11 -m1 0.1 -M1 50.0 -A1 0.01 -d1 200 -iP 0 -en 0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 1.5 -H1 0.0 -H6 5.0 -H9 10.0 -H10 2.0    
   
   . ~$ ./DIFFUSION_MR -y0 7 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -I0 16 -m0 0.1 -M0 100.0 -A0 0.01 -d0 200  -I1 17 -m1 0.01 -M1 10.0 -A1 0.01 -d1 200 -iP 0 -en 0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 1.5 -H1 1.0 -H6 1.0 -H9 10.0 -H10 2.0

   Log Scale Representation: x axis log scale (Nu) , y axis linear scale (alpha) 
   . ~$ ./DIFFUSION_MR   -y0 7 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -I1 16 -m1 1.0 -M1 10.0 -A1 0.01 -d1 200  -I0 17 -m0 0.01 -M0 10.0 -A0 0.01 -d0 200 -iP 0 -en 0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 2.5 -H1 1.0 -H6 1.0 -H9 10.0 -H10 2.0

   . ~S ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -G0 1 -G1 1 -sT 1.0E-06 -sN 300 -sP 2 -I1 16 -m1 1.0 -M1 10.0 -A1 0.01 -d1 200  -I0 17 -m0 0.01 -M0 10.0 -A0 0.01 -d0 200 -iP 0 -en 0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 2.5 -H1 1.0 -H6 1.0 -H9 10.0 -H10 2.0 -H11 0.0 -H12 0.0
*/

gsl_rng * r; /* Global generator defined in main.c */

int main(int argc, char **argv)
{
  int i, k, key;
  char * pF; 
  Parameter_Table Table;
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

  Parameter_Model * Par_Values = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (&Table, Par_Values);
  printf(" Parameter_Model structure 'Par_Values' has been correctly allocated and initiated\n");
  
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

  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( &Table, Par_Values) ;


  /* This vector is used by Fixed_Points_All(...) function, 
     which is called by Function_to_Type_of_Stability 
  */
  Table.Vector_Model_Variables_Stationarity = (double *)calloc(Table.MODEL_STATE_VARIABLES,
							       sizeof(double) );
  /* B E G I N : Main Function Call -------------------------------------------------*/  
  double * W_GRID = (double *)malloc( No_of_POINTS_1 * No_of_POINTS_2 * sizeof(double) );
  int Status =  generic_Function_Parameter_2Dim_Scan(&Table, 
						     No_of_POINTS_1, Input_Parameter_1,
						     No_of_POINTS_2, Input_Parameter_2,
						     Function_to_Type_of_Stability, 
						     W_GRID, "Coexistence_Condition.dat");
  /*   E N D : ----------------------------------------------------------------------*/
  free(Table.Vector_Model_Variables_Stationarity);
  
#if defined CPGPLOT_REPRESENTATION
  /* BEGIN : 2D GRID cpgplot representation */
      /*********************************************************************/
      Table.CPG->X_label   = Table.Symbol_Parameters_Greek_CPGPLOT[Input_Parameter_1]; 
      Table.CPG->Y_label   = Table.Symbol_Parameters_Greek_CPGPLOT[Input_Parameter_2]; 
      /*********************************************************************/
      
      Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_min );
      Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_MAX );

      /* X Axis: Log Scale */
      Table.CPG->ORIGIN_X    = log10(Value_0);
      Table.CPG->X_Dimension = log10(Value_1) - log10(Value_0);
      /* Table.CPG->ORIGIN_X    = Value_0;            */
      /* Table.CPG->X_Dimension = Value_1 - Value_0;  */
      
      Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_min );
      Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_MAX );

      /* Y Axis: Linear Scale */
      /* Table.CPG->ORIGIN_Y = log10(Value_0);                     */
      /* Table.CPG->Y_Dimension = log10(Value_1) - log10(Value_0); */
      Table.CPG->ORIGIN_Y = Value_0;
      Table.CPG->Y_Dimension = Value_1 - Value_0;
      
      Table.CPG->x_GRID  = No_of_POINTS_1; 
      Table.CPG->y_GRID  = No_of_POINTS_2;
      
      int Output_Variable  = Table.OUTPUT_VARIABLE_INDEX[0];
      Table.CPG->W_label   = Table.Output_Variable_Name[Output_Variable];

      Table.CPG->Title[0]='\0';  
      // pF = strcat(Table.CPG->Title, "Type of Stability Regimes (axes in log scale)");
      pF = strcat(Table.CPG->Title, "Type of Stability Regimes (X axis in log scale)");

      double * Color_States = (double *)calloc(4, sizeof(double));
      Color_States[0] = 3.0;  /* green  */ /* Only Resources                      */
      Color_States[1] = 2.0;  /* Red    */ /* Coexistence: non-damped oscillions  */
      Color_States[2] = 5.0;  /* Blue   */ /* Coexistence: damped oscillations    */
      Color_States[3] = 8.0;  /* Orange */ /* Coexistence: limit cycles           */
      
      int FIRST_PLOT = 0;
      double i_PLOT  = 0.0;
      C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___R_E_C_T_A_N_G_L_E_S( Table.CPG,
								    W_GRID, 
								    FIRST_PLOT,
								    4, Color_States,
								    0 );
      FIRST_PLOT = 1;
      C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___R_E_C_T_A_N_G_L_E_S___F_R_A_M_E(Table.CPG,
      									     W_GRID,
      									     FIRST_PLOT,
      									     4, Color_States,
      									     0 );
      free(Color_States);
#endif    
  /*   END : 2D GRID cpgplot representation */
  free(W_GRID);
  free(Par_Values);
  
  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */ 
#if defined CPGPLOT_REPRESENTATION
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif  

  #include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );
  free(Initial_Condition_Space);  
  free(Error_Space);  
 
  #include <include.Output_Variables.default.free.c>

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */
  
  printf("\nEnd of progam\n");
  return (0);
}
