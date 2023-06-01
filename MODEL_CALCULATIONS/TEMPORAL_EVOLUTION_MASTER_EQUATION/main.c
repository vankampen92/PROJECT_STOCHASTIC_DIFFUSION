/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2022 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <Include/MODEL.h>

#include "global.h"

gsl_rng * r; /* Global generator defined in main.c */

/* This code calculates the temporal evolution of the master equation corresponding to 
   several resource-consumer models 

   Compilation (see first makefile variable MODEL to specify a given model to compile):

   . ~$ make 

   Exectution:
   
   Single patch (-HM 1 -HX 2 -HY 1), and 3 species ---A, RA (and ARA, but only two dynamic
   variables, A and RA). Notice -H11 [Chi] -H12 [Eta]. If these two parameters are zero, there 
   is no triplet formation, and the feeding model is HOLLING Type II (-H9 [Alpha] -H10 [Nu]). 
   
   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 20 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 14 -HK 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H9 8.0 -H10 2.0 -H11 50.0 -H12 0.5 -Hp1 0.4 -Hp2 0.5 -HN 20 -tE 0.1

   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 20 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 14 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H9 2.5 -H10 10.0 -H11 100.0 -H12 1.0 -Hp1 0.3725 -Hp2 0.5 -HN 20 -tE 0.1 
   
   . ~$ ./DIFFUSION_HII_1D -y0 12 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 1 -v0 0 -G0 1 -G1 1 -tn 50 -t0 0.0 -t1 1.5 -t4 0 -tR 10000 -xn 0 -xN 20.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 14 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H9 2.5 -H10 1.0 -H11 0.0 -H12 0.0 -Hp1 0.3750 -Hp2 0.5 -HN 20

   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 10 -t0 0.0 -t1 15 -t4 0 -tR 10000 -xn 0 -xN 40.0 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 40 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H9 2.5 -H10 10.0 -H11 100.0 -H12 1.0 -Hp1 0.3 -Hp2 0.05 -HN 40 -tE 0.1

   -HuR -HuC are the jumping rates
   -H0  -H5  are the external immigration (Lambda_R_0 and Lambda_C_0)
   -H9  and -H10 are the Alpha_C_0 and Nu_C_0  Holling Type II model parameters.
   -H11 and -H12 are the Chi_C_0 and Eta_C_0 the Beddington_DeAngelis parameters.  
   
  -Hp1: Resource Carrying Capacity Fraction   
  -Hp2: No of Free Predator a Time 0 Fraction 
  -HN: No_of_INDIVIDUALS (TOTAL No of CONSUMERS)

   See denition_OutPut_Variables.c to understand the difference between Genuine Output Variables
   and plain model variables.

   More examples in ./command_line_examples.txt.
*/

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
  Table.CPG     = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  Table.CPG_STO = A_C_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T ( 1, SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME );
#endif
  
  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition:  
   */
  /* This definition is contingent to TYPE of MODEL at work from the pre-defined family of models:
     DIFFUSION_BD_2D, DIFFUSION_HII_1D, ... 
  */
  /* double p_1;         */ /* -Hp1 */ /* Resource Carrying Capacity Fraction */ 
  /* double p_2;         */ /* -Hp2 */ /* See below the definition of the
                                       /* TOTAL_No_of_FREE_CONSUMERS_TIME_0 */
  
  // void Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
  
  Table.TOTAL_No_of_RESOURCES  = (int)(Table.p_1 * (double)Table.K_R);
  Table.TOTAL_No_of_CONSUMERS  = Table.No_of_INDIVIDUALS;  /* -HN 20 or -HN 40 as an input argument */ 

  assert(Table.p_2 <= 1.0 && Table.p_2 >= 0.0);  // 
  assert(Table.p_1 <= 1.0 && Table.p_1 >= 0.0);  // Fractions!!!  
  
  Table.TOTAL_No_of_FREE_CONSUMERS_TIME_0     = (int)(Table.p_2*(double)Table.TOTAL_No_of_CONSUMERS);
  Table.TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = Table.TOTAL_No_of_CONSUMERS - Table.TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  Table.TOTAL_No_of_FREE_CONSUMERS            = (int)(Table.p_2*(double)Table.TOTAL_No_of_CONSUMERS);
  Table.TOTAL_No_of_HANDLING_CONSUMERS        = Table.TOTAL_No_of_CONSUMERS - Table.TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  /* END ----------------------------------------------------------------------------
     This initial Condition involves no triplets at time t = 0.0 because the sum of states
     should add up the TOTAL No of CONSUMERS 
  */
  /* Deterministic Time Dynamics */
  Parameter_Values_into_Parameter_Table(&Table);   /* This is to make sure the same
						      parameter set as defined through
						      either the command line or the 
						      default files is used!!! 
						   */
  M_O_D_E_L( &Table );

#if defined STOCHASTIC_REALIZATIONS
  M_O_D_E_L___S_T_O( &Table );
#endif

  // Some models does no have a stochastic master equation
  // counter-part implemented yet! 
  // At the moment, only DIFFUSION_BD_2D and DIFFUSION_HII_1D do
#ifdef DIFFUSION_BD_2D
  /* Stochastic Master Equation Time Evolution */
  Parameter_Values_into_Parameter_Table(&Table);   /* This is to make sure the same
						                                          parameter set as defined through
						                                          either the command line or the 
						                                          default files is used!!! 
						                                       */
  M_O_D_E_L___M_E( &Table );
#elif defined DIFFUSION_HII_1D  
  /* Stochastic Master Equation Time Evolution */
  Parameter_Values_into_Parameter_Table(&Table);   /* This is to make sure the same
						                                          parameter set as defined through
						                                          either the command line or the 
						                                          default files is used!!! 
						                                       */
  M_O_D_E_L___M_E( &Table );
#elif defined DIFFUSION_HII_nD  
  /* Stochastic Master Equation Time Evolution */
  Parameter_Values_into_Parameter_Table(&Table);   /* This is to make sure the same
						                                          parameter set as defined through
						                                          either the command line or the 
						                                          default files is used!!! 
						                                       */
  M_O_D_E_L___M_E( &Table ); 
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

  /* BEGIN : De-allocating All Memmory * * * * * * * * * * * */
#include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );

#include <include.Initial_Conditions.default.free.c>

#include <include.Output_Variables.default.free.c>

#include <include.Time_Dependence_Control.default.free.c>
  if (TYPE_of_TIME_DEPENDENCE == 0) T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( &Time, &Table );
  else                        Time_Dependence_Control_Free( &Time_Dependence, &Table );

#if defined CPGPLOT_REPRESENTATION
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG_STO, SUB_OUTPUT_VARIABLES );
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */

  printf("\nEnd of progam\n");
  return (0);
}
