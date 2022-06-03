/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2018 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <Include/MODEL.h>

#include "global.h"

gsl_rng * r; /* Global generator defined in main.c */

/* This code calculates the stochastic and determinisitic temporal evolution of several resource-consumer models 

   Compilation (see makefile variable MODEL). Some compilation commands as example:

   . ~$ make
   . ~$ make STATIONARITY=STATIONARY_POINT_REPRESENTATION MODEL=DIFFUSION_STOLLENBERG_4D
   . ~$ make STATIONARITY=NON_STATIONARY_POINT_REPRESENTATION MODEL=DIFFUSION_STOLLENBERG_4D

   Exectution:
   
   1 species examples:                                       (OUTPUT_VARIABLES_GENUINE will be 4) 
   . ~$ ./DIFFUSION -y0 0 -y2 1 -HS 1 -HM 100 -HX 10 -HY 10 -Hu 0.1 -n 2 -v0 0 -v1 60 -G0 1 -G1 2 -tn 50 -t0 0.0 -t1 50.0 -t4 0 -tR 10 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 50.0 -G5 1 -G6 0.0 -G7 6000.0

   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 -Hu 0.1 -n 1 -v0 5054 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 200.0 -t4 0 -tR 10 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 50.0 -G5 1 -G6 0.0 -G7 1100.0

   2 species examples:                                       (OUTPUT_VARIABLES_GENUINE will be 5) 
   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 2 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -n 1 -v0 10105 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 30.0 -t4 0 -tR 5 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 30.0 -G5 1 -G6 0.0 -G7 1100.0

   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 2 -HM 100 -HX 10 -HY 10 -Hu 0.5 -n 1 -v0 115 -G0 1 -G1 1 -tn 20 -t0 0.0 -t1 5.0 -t4 0 -tR 5 -xn 0 -xN 98 -HN 98 -G2 1 -G3 0.0 -G4 5.0 -G5 1 -G6 0.0 -G7 100.0

   3 species example:                                        (OUTPUT_VARIABLES_GENUINE will be 6) 
   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 3 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -n 1 -v0 15156 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 30.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 30.0 -G5 1 -G6 0.0 -G7 1100.0

   1 species (with external immigration and death) example:  (OUTPUT_VARIABLES_GENUINE will be 4) 
   .~$ ../DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.1 -n 1 -v0 5054 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0

   2 species (with external immigration and death) example:  (OUTPUT_VARIABLES_GENUINE will be 5) 
   .~$ ./DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 2 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.1 -n 1 -v0 10105 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0
   Notice that the 2nd species (in green) does not die or externally immigrate.

   .~$ ./DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 2 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.5 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0
   Notice that the 2nd species (in green) does not die or externally immigrate.

   10 species (with external immigration and death) example (OUTPUT_VARIABLES_GENUINE will be 13) 
   .~$ ./DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 10 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.5 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0

   Important notice: MODEL.h contains a 'define' of No_of_RESOURCES_MAXIMUM. If you overcome that
   limit, program crashes. 

   4 species (4 different species ---R, A, RA, ARA---, therefore OUTPUT_VARIABLES_GENUINE is 7) 
   .~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100  -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0

   MODEL = DIFFUSION_1R1C
   Single patch (-HM 1 -HX 1 -HY 1), and 3 species ---R, A, RA. Notice -H11 [Chi_C_0] -H12 [Eta_C_0]. If these two parameters are zero, no triplet formation, and the dynamics is equivalent to a 3D system, with only three local model variables.
   . ~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 -tn 100 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -H6 0.5 -H9 10.0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 2.5 -H1 1.0 -H6 1.0 -H9 8.0 -H10 2.0 -H11 0.0 -H12 0.0

   MODEL = DIFFUSION_STOLLENBERG_3D
   Single patch (-HM 1 -HX 1 -HY 1), and 3 species ---R, A, RA. The dynamics do not include triplet formation. It is a 3D system, with only three local model variables (R, A, RA).
   . ~$ ./DIFFUSION_STOLLENBERG_3D -y0 10 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 -tn 100 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -HK 2000 -H1 1.0 -H6 0.5 -H9 10.0 -H10 2.0 -H4 3.5 -H17 1.0 
   See denition_OutPut_Variables.c to understand the difference between Genuine Output Variables
   and plain model variables.

   MODEL = DIFFUSION_STOLLENBERG_4D
   Single patch (-HM 1 -HX 1 -HY 1), and 4 species ---RP, R, A, RA. The dynamics do not include triplet formation. It is a 4D system, with only four local model variables (RP, R, A, RA). 
   . ~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 2000 -H20 20.0 -H1 1.0 -H3 5.0 -H6 0.5 -H9 20.0 -H10 2.0 -H4 5.0 -H17 1.0 
   This example produces a limit cycle. 
  
   . ~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 -tn 200 -t0 0.0 -t1 75.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 75.1 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 2000 -H20 20.0 -H1 1.0 -H3 5.0 -H6 0.5 -H9 5.0 -H10 2.0 -H4 5.0 -H17 1.0
   This example produces damped oscillations

   See denition_OutPut_Variables.c to understand the difference between Genuine Output Variables
   and plain model variables.
   
   -HuR -HuC are the jumping rates
   -H0 -H2 -H5  are the external immigration (Lambda_R_0, Lambda_R_1 and Lambda_C_0)
   -H20 is the establishment rate 
   -H1  -H3  -H6 are the death rates (Delta_R_0, Delta_R_1 for propagules, and Delta_C_0 for both searching and handling consumers)
   -H9  and -H10 are the Alpha_C_0 and Nu_C_0  Holling Type II model parameters 
   -H4  and -H17 are the production rates of propagules (Beta_R) and searching animals (Beta_C), respectively.  
    
   MacArthur and Rosenzweig (two species 3D, R, A, RA):
   .~$ ./DIFFUSION_MR -y0 7 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 -tn 100 -t0 0.0 -t1 40.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 -G2 1 -G3 0.0 -G4 40.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -HK 2000.0  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 25.0 -H1 0.0 -H6 5.0 -H9 17.0 -H10 5.0

   Coarse-grained 2D system (two species 1R and 1C):
   .~$ ./DIFFUSION_1R1C_2D -y0 4 -y2 1 -HS 1 -HM 4 -HX 2 -HY 2 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 50 -HN 50 -G2 1 -G3 0.0 -G4 40.0 -G5 1 -G6 0.0 -G7 100.0 -H0 0.01 -H5 0.01 -H6 0.5

   MODEL = DIFFUSION_BD_2D
   Feeding experiments at a constant number of total consumers: 
   .~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 50 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 20 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H9 2.5 -H10 10.0 -H11 100.0 -H12 1.0 -Hp1 0.3725 -Hp2 0.5 -HN 20

   More examples in ./command_line_examples.txt.
*/

void Common_Initial_Condition_Command_Line_Arguments_into_Table(Parameter_Table *Table);

int main(int argc, char **argv)
{
  int i;
  Parameter_Table Table;
  Time_Control Time;
  Time_Dependence_Control Time_Dependence;
  P_ARG = &Table;

#include "default.c"

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

#include <include.Output_Variables.default.aux.c>
  
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C(   &Table );
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( &Table, Index_Output_Variables );
  printf(" Parameter_Table structure has been correctly allocated and initiated\n");


  /* B E G I N : Reserving memmory for Parameter Space */
#include <include.Parameter_Space.default.aux.c>
  if( No_of_PARAMETERS == Table.TOTAL_No_of_MODEL_PARAMETERS ) {
    /* Full parameter space is in place. See also Model_Variables_Code.c */
    for(i=0; i<Table.TOTAL_No_of_MODEL_PARAMETERS; i++) Index[i] = Table.Index[i];
    No_of_PARAMETERS = Table.TOTAL_No_of_MODEL_PARAMETERS;
  }
  Parameter_Space * Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  Parameter_Space_Alloc( Space, No_of_PARAMETERS, d);
  Parameter_Space_Initialization( Space, No_of_PARAMETERS, TOLERANCE, MAX_No_of_ITERATIONS,
    d, Index, Ranges, Acc);
  Table.S = Space;
  printf(" Parameter_Space structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */

#include <gsl_random_number_Setup.c>
  // #if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); Press_Key();
  /*   END: Checking Random Number Generator Setup */
  // #endif

  if (TYPE_of_TIME_DEPENDENCE == 0) {
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

    int No_of_EMPIRICAL_TIMES = I_Time;
    // Number of columns in the data files of time-dependent parameters
    Time_Dependence_Control_Upload(&Time, &Time_Dependence, &Table,
				   I_Time, No_of_EMPIRICAL_TIMES,
				   TIME_DEPENDENT_PARAMETERS, TYPE_of_TIME_DEPENDENCE,
				   TYPE_0_PARAMETERS, TYPE_1_PARAMETERS, TYPE_2_PARAMETERS,
				   No_of_COVARIATES,
				   dependent_parameter, forcing_pattern,
				   "File_of_Covariates.dat", Name_of_FILE[0] );
    printf(" Both Time_Control and Time_Dependence_Control structures have been\n");
    printf(" correctly allocated and set up\n");
  }

#if defined CPGPLOT_REPRESENTATION
  Table.CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  Table.CPG_STO = A_C_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T (0,
							 SUB_OUTPUT_VARIABLES, I_Time,
							 0, CPG_DRIVER_NAME);
  printf(" Two Parameterh_CPGPLOT plotting structures have been correctly allocated and initiated\n");
  printf(" These will open two windows (or two ploting devices of the same kind)\n");
  printf(" Table.CPG will store deterministic dynamic variables to plot\n");
  printf(" Table.CPG_STO will store stochastic dynamic variables to plot\n");
  printf(" As a consquence, deterministic and stochastic dynamics can be plotted\n");
  printf(" on the same device to compare (as it is done here, indicated by the first\n");
  printf(" input argument (0) of the A_Ch_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T function).\n");
  printf(" Alternatively, two different devices (two different pdf files, for instance)\n");
  printf(" can be used, if required (1).\n");
#endif

  if(Table.TYPE_of_MODEL == 12 || Table.TYPE_of_MODEL == 13 || Table.TYPE_of_MODEL == 14) 
    /* Models where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
    
  /* Deterministic Time Dynamics */
  Parameter_Values_into_Parameter_Table(&Table);
  M_O_D_E_L( &Table );
  
  // Some models (such as DIFFUSION_1R1C_2D) do no have a stochastic
  // counter-part implemented yet!
#ifndef DIFFUSION_1R1C_2D
#ifndef DIFFUSION_DRAG
#ifndef DIFFUSION_VRG
#ifndef DIFFUSION_MR
  /* Stochastic Time Dynamics: A number of stochastic realizations     */
  Parameter_Values_into_Parameter_Table(&Table);
  M_O_D_E_L___S_T_O( &Table );
#endif
#endif
#endif
#endif
  /* BEGIN : -------------------------------------------------------------------------
   */
  char boundary_File[80];
  sprintf(boundary_File, "boundary_Model_Parameter.c");
  write_Parameter_Table___RANGES___VALUES___LATEX ( "Latex_Parameter_Table.tex",
                                                    boundary_File,
                                                    &Table,
                                                    Space->P_MAX->data,
                                                    Space->P_min->data, Space->No_of_PARAMETERS);
  /*  END : ------------------------------------------------------------------------*/

  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */
#if defined CPGPLOT_REPRESENTATION
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG_STO, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif

#include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );

#include <include.Initial_Conditions.default.free.c>

#include <include.Output_Variables.default.free.c>

#include <include.Time_Dependence_Control.default.free.c>
  if (TYPE_of_TIME_DEPENDENCE == 0) T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( &Time, &Table );
  else                        Time_Dependence_Control_Free( &Time_Dependence, &Table );

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */

  printf("\nEnd of progam\n");
  return (0);
}

void Common_Initial_Condition_Command_Line_Arguments_into_Table(Parameter_Table *Table)
{
      /* BEGIN : -------------------------------------------------------------------------
       * Definition Initial Condition:  
       */
      /* This definition is contingent to TYPE of MODEL at work from the pre-defined family 
	 of models:
	 DIFFUSION_BD_2D, DIFFUSION_HII_1D, ... 
      */
      /* double p_1;         */ /* -Hp1 */ /* Resource Carrying Capacity Fraction */ 
      /* double p_2;         */ /* -Hp2 */ /* See below the definition of the     */
                                           /* TOTAL_No_of_FREE_CONSUMERS_TIME_0   */
      Table->TOTAL_No_of_RESOURCES  = (int)(Table->p_1 * (double)Table->K_R);
      Table->TOTAL_No_of_CONSUMERS  = Table->No_of_INDIVIDUALS;  /* -HN 20 as input argument */ 

      assert(Table->p_2 <= 1.0 && Table->p_2 >= 0.0);  // 
      assert(Table->p_1 <= 1.0 && Table->p_1 >= 0.0);  // Fractions!!!  
  
      Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = (int)(Table->p_2*(double)Table->TOTAL_No_of_CONSUMERS);
      Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = Table->TOTAL_No_of_CONSUMERS - Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
      /* END ----------------------------------------------------------------------------
	 This initial Condition involves no triplets at time t = 0.0 because the sum of states
	 should add up the TOTAL No of CONSUMERS 
      */
}
  
