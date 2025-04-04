#include "HEADERS.h"
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                            David Alonso, 2022 (c)                         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define HUMAN_CONSTANT_POPULATION

#define INTEGER_CODE_FOR_TIME_DIMENSION 0

/* Model parameters are drawn from a pool. Only the most complex model defined might use them 
   all. Simpler models only use a subset from the whole pool. The maximum number of model 
   parameters for each specific model is defined in the corresponding file: 

   MODEL_DEFINE_MAX_VARIABLES_[MODEL].h, 

   which defines the maxim dinension of the whole parameter space for that particular model 
*/
#define MAX_No_of_CONFIGURATIONAL_STATES 100000000 /* Max No of Eqs in the Master Equation */   

#define No_of_TDC_FUNC_AUX_PARAM_MAX 10 /* Maximum No of Time Dependence Function Auxiliary 
                                           Parameters required to define the potential 
                                           functional time dependence for each model parameter 
                                        */
#define DEPENDENT_PARAMETERS_MAXIMUM 30 /* Maximum number of potentially time-dependent 
                                           forcing parameters */
#define MODEL_PARAMETERS_MAXIMUM 30     /* Maximum No of MODEL (input) PARAMETERS */
                                        /* The total number of parameters in all 
					                                 model-paramater-related assign functions.
					                                 This is the whole parameter pool from which a 
					                                 parameter subspace can be defined for 
					                                 optimization searches, parameter scans,   
					                                 and each specific model
					                              */
#define ME_n_DIMENSION_MAXIMUM 20       /* Maximum dimension of the probability distribution 
                                        // in a master equation */
#ifdef DIFFUSION_1RnC_E
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_1RnC_E.h"
#elif defined DIFFUSION_1R1C
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_1R1C.h"
#elif defined DIFFUSION_MR
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_MR.h"
#elif defined DIFFUSION_S_RESOURCES
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_S_RESOURCES.h"
#elif defined DIFFUSION_1R1C_2D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_1R1C_2D.h"
#elif defined DIFFUSION_1R1C_2D_STO_4D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_1R1C_2D_STO_4D.h"
#elif defined DIFFUSION_DRAG
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_DRAG.h"
#elif defined DIFFUSION_VRG
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_VRG.h"
#elif defined DIFFUSION_HII_2D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_HII_2D.h"
#elif defined DIFFUSION_HII_nD
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_HII_nD.h"
#elif defined DIFFUSION_STOLLENBERG_3D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_STOLLENBERG_3D.h"
#elif defined DIFFUSION_HII_AC_2D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_HII_AC_2D.h"
#elif defined DIFFUSION_HII_1D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_HII_1D.h"
#elif defined DIFFUSION_BD_2D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_BD_2D.h"
#elif defined DIFFUSION_BD_3D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_BD_3D.h"
#elif defined DIFFUSION_STOLLENBERG_4D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_STOLLENBERG_4D.h"
#elif defined DIFFUSION_AZTECA_4D
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_AZTECA_4D.h"
#elif defined DIFFUSION_AZTECA_4D_0
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_AZTECA_4D_0.h"
#elif defined DIFFUSION_AZTECA_4D_1
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_AZTECA_4D_1.h"
#elif defined DIFFUSION_ECOEVO_PLANTS
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_ECOEVO_PLANTS.h"
#elif defined DIFFUSION_ECO_PLASMIDS
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_ECO_PLASMIDS.h"
#elif defined DIFFUSION_ECO_1B1P
#include "MODEL_DEFINE_MAX_VARIABLES_DIFFUSION_ECO_1B1P.h"  
#endif

typedef struct totalRateinfo
{
  double Total_Rate;
  double max_Probability;
  double Stochastic_Time;
  double Reusable_Random_Number; 
}Stochastic_Rate;

#ifndef BINARY_TREE_STRUCTURE 
  #define BINARY_TREE_STRUCTURE
  typedef struct treenode
  {
    /* In this implementation, only true tree leaves are ordered      */
    /* The order of internal nodes is not initialized                 */
    int    level; /* Tree Level (0: root; n: leaves)                  */
    int    order; /* Order within the leaf level, from 0 to 2^n - 1   */
    int    index; /* Index (when using the binary tree as a priority queu) */

    double value;

    struct treenode * left;
    struct treenode * right;
    struct treenode * parent;
  }treenode;
#endif
/* More Stucture Defintions:                               */
#include "MODEL_Trend_Control_STRUCT_DEF.h"
#include "MODEL_Time_Control_STRUCT_DEF.h"
#include "MODEL_Time_Dependence_Control_STRUCT_DEF.h"
#include "MODEL_Parameter_Model_STRUCT_DEF.h"
#include "MODEL_Strain_STRUCT_DEF.h"
#include "MODEL_Plasmid_STRUCT_DEF.h"
#include "MODEL_Community_STRUCT_DEF.h"
#include "MODEL_Parameter_Space_STRUCT_DEF.h"
#include "MODEL_Master_Equation_STRUCT_DEF.h"
#include "MODEL_Parameter_Table_STRUCT_DEF.h"
#include "MODEL_Generic_Root_Data_STRUCT_DEF.h"
#include "MODEL_Observed_Data_STRUCT_DEF.h"
#include "MODEL_Parameter_Fitting_STRUCT_DEF.h"
/*  * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <Time_Control.h>
#include <Observed_Data.h>
#include <Parameter_Model.h>
#include <Parameter_Table.h>
#include <Community.h>
#include <In_Out_Migration_Functions.h>

/* Auxliary Libraries at ./Library common directory */
#include <GSL_Optimization.h>
#include <IO_Procedures.h>
#include <GSL_stat.h>

/* ...and more prototype definitions */
#include <Definition_Error_Model/Error_Library.h>
#include <Definition_Fixed_Points/Fixed_Points_All.h>
#include <Definition_Master_Equation/Master_Equation_Functions.h>
#if defined DIFFUSION_HII_nD
  #include <Definition_Master_Equation/DIFFUSION_HII_nD/Model_Parameters_Master_Equation.h>
#endif 

/* More Auxiliary Functions Prototype defintions */
#include "main.H"
#include <assign.h>
#include <Model_Variables_Code.h>
#include <definition_OutPut_Variables.h>
#include <definition_OutPut_Variables_Functions.h>

// #include <R_0.h>
#include <fprintf_Model_Parameters.h>
#include <fprintf_Output_Variables.h>
#include <Latex_Write_Parameter_Table.h>

// #include <Fixed_Points/parameter_Simple_Scan.h>
// #include <Fixed_Points/bifurcation_Diagram.h>
// #include <Fixed_Points_Functions/Fixed_Points_Functions.h>

// #include <generic_Function_Parameter_Scan.h>
// #include <generic_Subregion_Parameter_Space.h>
// #include <generic_Random_Parameter_Space.h>
// #include <generic_Random_Parameter_Space_Plus.h>

// #include <eigen_Functions.h>
// #include <Eigenvalue_Calculation.h>

#include <Definition_Time_Trend_Dependence/Time_Dependence_Control.h>
#include <Definition_Time_Trend_Dependence/Trend_Control.h>
#include <Definition_Parameter_Space/Parameter_Space.h>
#include <Definition_Fitting_Structure/Fitting_Structure.h>
#include <Definition_Stochastic_Realizations/Stochastic_Time_Dynamics.h>
#include <Definition_Numerical_Integration/deterministic_time_dynamics.h>
#include <Definition_Numerical_Integration/numerical_Integration_Driver.h>
#include <Definition_Numerical_Integration/Initial_Conditions_Numerical_Integration.h>
#include <Definition_Numerical_Integration/ODE_Definitions/ODE_Definitions.h>
#include <Definition_Master_Equation/Master_Equation_Functions.h>

#if defined CPGPLOT_REPRESENTATION
/* Header file for Parameter Table dependent CPGPLOT plotting auxiliary functions */
#include <CPGPLOT_Parameter_Table/CPGPLOT___X_Y___Parameter_Table.h>
#include <CPGPLOT_Parameter_Table/CPGPLOT___GRID___Parameter_Table.h>
#else
#define MAX_No_of_FILES 20
#endif
