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

   . ~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                          -n 2 -v0 0 -v1 1 -tn 10 -t0 0.0 -t1 15 -t4 0 -tR 10000 -tE 0.1 \
                          -xn 0 -xN 40.0 \
                          -G0 1 -G1 2 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 40 \
                          -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 \
                          -HK 10000 -H9 2.5 -H10 10.0 -H11 100.0 -H12 1.0 \
                          -Hp1 0.3 -Hp2 0.05 -HN 40 

   -HuR -HuC are the jumping rates
   -H0  -H5  are the external immigration (Lambda_R_0 and Lambda_C_0)
   -H9  and -H10 are the Alpha_C_0 and Nu_C_0  Holling Type II model parameters.
   -H11 and -H12 are the Chi_C_0 and Eta_C_0 the Beddington_DeAngelis parameters.  
   
  -Hp1: Resource Carrying Capacity Fraction   
  -Hp2: No of Free Predator a Time 0 Fraction 
  -HN: No_of_INDIVIDUALS (TOTAL No of CONSUMERS)

  MODEL = DIFFUSION_HII_nD
   Feeding experiments at a constant number of total consumers (-HN 20): 
   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 
   			                  -n 3 -v0 0 -v1 1 -v2 2 \
			                    -G0 1 -G1 3 -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 20 \
			                    -tn 50 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 \
			                    -HuR 0.0 -HuC 0.0 -H5 0.0 \
			                    -HK 1000 -H0 5.0 -H2 1.0 -H9 2.5 -H10 10.0 \
			                    -Hp1 0.3725 -Hp2 0.5 -HN 20

   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \
                          -n 3 -v0 0 -v1 1 -v2 2 \
                          -G0 1 -G1 3 -G2 1 -G3 0.0 -G4 2.5 -G5 1 -G6 0.0 -G7 8 \
                          -tn 10 -t0 0.0 -t1 2.5 -t4 0 -tR 100 -xn 0 -xN 20.0 -tE 0.2 \
                          -HuR 0.0 -HuC 0.0 -H5 0.0 \
                          -HK 10000 -H0 5.0 -H2 1.0 -H9 10.5 -H10 10.0 \
                          -Hp1 0.3725 -Hp2 1.0 -HN 20 

   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 5 -HM 1 -HX 1 -HY 1 \
                          -n 5 -v0 0 -v1 1 -v2 2 -v3 3 -v4 4\
                          -G0 1 -G1 5 -G2 1 -G3 0.0 -G4 3.5 -G5 1 -G6 0.0 -G7 8 \
                          -tn 5 -t0 0.0 -t1 3.5 -t4 0 -tR 10000 -xn 0 -xN 20.0 -tE 0.2 \
                          -HuR 0.0 -HuC 0.0 -H5 0.0 \
                          -HK 10000 -H0 5.0 -H2 1.0 -H9 10.5 -H10 0.1 \
                          -Hp1 0.3725 -Hp2 1.0 -HN 20 
   
   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 10 -HM 1 -HX 1 -HY 1 \
                          -n 10 -v0 0 -v1 1 -v2 2 -v3 3 -v4 4 -v5 5 -v6 6 -v7 6 -v8 8 -v9 9\
                          -G0 2 -G1 5 -G2 1 -G3 0.0 -G4 3.5 -G5 1 -G6 0.0 -G7 8 \
                          -tn 10 -t0 0.0 -t1 3.5 -t4 0 -tR 10000 -xn 0 -xN 20.0 -tE 0.2 \
                          -HuR 0.0 -HuC 0.0 -H5 0.0 \
                          -HK 10000 -H0 5.0 -H2 1.0 -H9 10.5 -H10 0.1 \
                          -Hp1 0.3725 -Hp2 1.0 -HN 20

   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 5 -HM 1 -HX 1 -HY 1 \
                          -n 5 -v0 0 -v1 1 -v2 2 -v3 3 -v4 4 \
                          -G0 1 -G1 5 -G2 1 -G3 0.0 -G4 6.0 -G5 1 -G6 0.0 -G7 8 \
                          -tn 5 -t0 0.0 -t1 6.0 -t4 0 -tR 10000 -xn 0 -xN 20.0 -tE 0.2 \
                          -HuR 0.0 -HuC 0.0 -H5 0.0 \
                          -HK 10000 -H0 5.0 -H2 1.0 -H9 0.5 -H10 0.1 \
                          -Hp1 0.3725 -Hp2 1.0 -HN 20

   Relevant input arguments for model DIFFUSION_HII_nD:

   -HK  [N or Total Carrying Capacity (for resources)]
   -Hp1 [approx y_R = f_i * p1 * N: Different resource levels at time 0.0, with f_i = 0.5, 0.3, 0.2, but this can be changed] 
   -Hp2 [Total No of Free Consumers at time 0.0 = p2 * No_of_CONSUMERS]
   -HN 20 [Total No of CONSUMERS]
   -H9 and -H10 [Alpha_C_0 and Nu_C_0 for 1st Resource Type. Alpha_C_0 is the same for all resource types, but this can be changed.  
   -H0 and -H2 are Lambda_R_0 and Lambda_R1, which are overloaded to create extra handling rates (Nu), for the 2nd and 3rd resource type. 
   Notice that for this model -H5 should be always zero (No immigration of consumers: No of CONSUMERS is constant). 
   -H11 and H12 [Xhi_C_0 and Eta_C_0 are not relevant here. They are only in models with predator interference]. 
   
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
  printf("\n"); Print_Press_Key(1,0,".");
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
    Time_Dependence_Control_Alloc(&Time, &Time_Dependence, &Table, I_Time, 
                                  TIME_DEPENDENT_PARAMETERS, No_of_COVARIATES);

    int No_of_EMPIRICAL_TIMES = I_Time;
    // Number of columns in the data files of time-dependent parameters
    Time_Dependence_Control_Upload( &Time, &Time_Dependence, &Table,
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
  if(Table.TYPE_of_MODEL == 12 || Table.TYPE_of_MODEL == 13 || Table.TYPE_of_MODEL == 14) 
    /* Models where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);

  if(Table.TYPE_of_MODEL == 16) { // DIFFUSION_HIIl_nD
    /* Also model where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    /* and they feed on multiple resources                      */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
    // Resetting_Alpha_Nu_Vectors (&Table);
    Resetting_Alpha_Nu_Vectors_Constant (&Table); /* Alpha and Nu parameter all the same */
    Resetting_Multiresource_Levels (&Table);      /* Creating Vector of Thetas */
    Writing_Alpha_Nu_Theta_Vectors(&Table);  

    Resetting_Multiresource_Levels (&Table);  
    Writing_Alpha_Nu_Theta_Vectors(&Table);  
  }
  /* Deterministic Time Dynamics */
  Parameter_Values_into_Parameter_Table(&Table);
  M_O_D_E_L( &Table );

#if defined STOCHASTIC_REALIZATIONS
  M_O_D_E_L___S_T_O( &Table );
#endif

  // Some models does no have a stochastic master equation
  // counter-part implemented yet! 
  // At the moment, only DIFFUSION_BD_2D, DIFFUSION_HII_1D, 
  // and DIFFUSION_HII_nD do

  /* Stochastic Master Equation Time Evolution */
  Parameter_Values_into_Parameter_Table(&Table);   /* This is to make sure the same
						                                          parameter set as defined through
						                                          either the command line or the 
						                                          default files is used!!! 
						                                       */
#ifdef DIFFUSION_BD_2D
  M_O_D_E_L___M_E( &Table );
#elif defined DIFFUSION_HII_1D  
  M_O_D_E_L___M_E( &Table );
#elif defined DIFFUSION_HII_nD  
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
