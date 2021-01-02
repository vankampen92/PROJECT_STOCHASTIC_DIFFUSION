#include "HEADERS.h"
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                            David Alonso, 2021 (c)                         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define HUMAN_CONSTANT_POPULATION

#define INTEGER_CODE_FOR_TIME_DIMENSION 0

#define DEPENDENT_PARAMETERS_MAXIMUM 11   /* Maximum number of potentially forced parameters */

#define MODEL_PARAMETERS_MAXIMUM 11       /* Maximum No of MODEL (input) PARAMETERS */

#define MODEL_PARAMETER_SPACE_MAXIMUM 11  /* Maximum Dimension for Parameter Space */

#define No_of_CELLS_MAXIMUM 10000

#define No_of_RESOURCES_MAXIMUM 3             /* S *//* Number of local states */

#define MODEL_STATE_VARIABLES_MAXIMUM 30000   /* M */ /* Dimension maximum: 2 * 10000 cells */
                                            /* S * M */
#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 6  /* Number Output Variables            */
                                            /* (other than MODEL_STATE_VARIABLES) */
                                            /* 3 (Resources Maximum) + 3          */
                                            /* See definition_OutPut_Variables.c  */

#define OUTPUT_VARIABLES_MAXIMUM 30006  /* MODEL_STATE_VARIABLES_MAXIMUM +
					       OUTPUT_VARIABLES_GENUINE_MAXIMUM */

typedef struct totalRateinfo
{
  double Total_Rate;
  double max_Probability;
  double Stochastic_Time;
}Stochastic_Rate;

#include "MODEL_Trend_Control_STRUCT_DEF.h"

#include "MODEL_Time_Control_STRUCT_DEF.h"

#include "MODEL_Time_Dependence_Control_STRUCT_DEF.h"

#include "MODEL_Parameter_Model_STRUCT_DEF.h"

#include "MODEL_Community_STRUCT_DEF.h"

#include "MODEL_Parameter_Space_STRUCT_DEF.h"

#include "MODEL_Parameter_Table_STRUCT_DEF.h"

#include "MODEL_Generic_Root_Data_STRUCT_DEF.h"

#include "MODEL_Observed_Data_STRUCT_DEF.h"

#include "MODEL_Parameter_Fitting_STRUCT_DEF.h"

#include <Time_Control.h>

#include <Parameter_Model.h>

#include <Parameter_Table.h>

#include <Community.h>

#include <In_Out_Migration_Functions.h>

/* Auxliary Libraries at ./Library common directory */
#include <GSL_Optimization.h>

#include <IO_Procedures.h>

#include <GSL_stat.h>

#include <Definition_Error_Model/Error_Library.h>

/* Auxiliary Functions */
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
#include <Definition_Stochastic_Realizations/Stochastic_Time_Dynamics.h>
#include <Definition_Numerical_Integration/deterministic_time_dynamics.h>
#include <Definition_Numerical_Integration/numerical_Integration_Driver.h>
#include <Definition_Numerical_Integration/Initial_Conditions_Numerical_Integration.h>
#include <Definition_Numerical_Integration/ODE_Definitions/ODE_Definitions.h>

#if defined CPGPLOT_REPRESENTATION
/* Header file for Parameter Table dependent CPGPLOT plotting auxiliary functions */
#include <CPGPLOT_Parameter_Table/CPGPLOT___X_Y___Parameter_Table.h>
#else
#define MAX_No_of_FILES 20
#endif
