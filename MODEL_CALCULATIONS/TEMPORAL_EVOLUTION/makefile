#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/*                           David Alonso, 2021 (c)                          */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#
# Note: If rather than using the recursive make at the project root directory, 
# you compile only from here, compilation may fail when calling make. 
# It means that you need to export these two environment variables:
#
# . ~$ export PROJECT_DIR=${HOME}/PROJECT_AB-COVID19_MODEL/
# . ~$ export CPGPLOT=${HOME}/CPGPLOT/
#
# Alternatively, you can also define these environmental variables here in
# the make file:
CPGPLOT=${HOME}/CPGPLOT/
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
INCLUDE_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/Include/
TIME_TREND_DIR=${PROJECT_DIR}Definition_Time_Trend_Dependence/
IO_DIR=${PROJECT_DIR}Library/IO_Library/
OPT_DIR= ${PROJECT_DIR}Library/Optimization_Library/
ERROR_DIR=${PROJECT_DIR}Definition_Error_Model/
FIXED_POINTS_DIR=${PROJECT_DIR}Definition_Fixed_Points/###########
WKD=${PROJECT_DIR}MODEL_CALCULATIONS/TEMPORAL_EVOLUTION/#
DET_DIR=${PROJECT_DIR}MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/#
C_NR=${PROJECT_DIR}Library/C_NumRecip/#############
GSL_DA=${PROJECT_DIR}Library/GSL_RANDOM/##########
HEAD.H =${PROJECT_DIR}Include/MODEL.h##############
# Mac users, please add this directore to the list 
# of library paths:
# MAC_OPT_DIR=/opt/local/include/
INCLUDE_DIR=${PROJECT_DIR}Include/#################
# Directories where libraries are located 
# the linker will look for at the final linking step, 
# after all object files, *.o, are successfully compiled,
# to create the exectuable code file 
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}###################
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}library/###########
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}library_basic_UTL/#
EXTRA_LIBRARY_PATH+=-L${C_NR}######################
EXTRA_LIBRARY_PATH+=-L${GSL_DA}####################
EXTRA_LIBRARY_PATH+=-L${IO_DIR}
EXTRA_LIBRARY_PATH+=-L${OPT_DIR}
EXTRA_LIBRARY_PATH+=-L${ERROR_DIR}
# Directories where headers files are located:
# (required at the precompilation step, when building 
# the files to compile and create the *.o objects with 
# the suitable header files)
# Mac user, please add:
# EXTRA_INCLUDE_PATH =-I${MAC_OPT_DIR}###############
# Linux users, simply: 
EXTRA_INCLUDE_PATH+=-I${INCLUDE_DIR}###############
EXTRA_INCLUDE_PATH+=-I${PROJECT_DIR}###############
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}###################
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library/###########
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library_basic_UTL/#
EXTRA_INCLUDE_PATH+=-I${C_NR}######################
EXTRA_INCLUDE_PATH+=-I${GSL_DA}####################
EXTRA_INCLUDE_PATH+=-I${IO_DIR}
EXTRA_INCLUDE_PATH+=-I${OPT_DIR}
EXTRA_INCLUDE_PATH+=-I${ERROR_DIR}
EXTRA_INCLUDE_PATH+=-I${FIXED_POINTS_DIR}
##############################################################################################################
## C O M P I L A T I O N   O P T I O N   C O N T R O L I N G   N U M E R I C A L   R E C I P E S
NR=NUMERICAL_RECIPES###############################
##############################################################################################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   S T O C H A S T I C   D Y N A M I C S
## STOCHASTIC_OPTIMIZATION requires de definition of the Events Adjancy Matrix telling which events
##                         influences the rate of occurrence of any other event, locally
## BINARY_TREE_OPTIMIZATION requires setting up a binary tree for all the spatial individual cells
## BINARY_TREE_SUPER_OPTIMIZATION requires setting up a bigger binary tree of every event at each cell 
STO=STOCHASTIC_OPTIMIZATION###################### Possibilities: 'STOCHASTIC_OPTIMIZATION' / 'NON'
BT=BINARY_TREE_SUPER_OPTIMIZATION#################### Possibilities: 
##################################################### (0) 'NON' (compulsary if STO=NON)
##################################################### (1) 'BINARY_TREE_OPTIMIZATION'
##################################################### (2) 'BINARY_TREE_SUPER_OPTIMIZATION'
##################################################### (3) 'PRIORITY_QUEU_SUPER_OPTIMIZATION'
RE=REUSE_RANDOM_NUMBER############################### (22) Possibilities;
#####################################################      'NON_REUSE_RANDOM_NUMBER'
#####################################################      'REUSE_RANDOM_NUMBER'
# The optimation option of reusing random numbers is only available in case 
# BT=BINARY_TREE_SUPER_OPTIMIZATION (otherwise does not make any sense)
##############################################################################################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   P R O G R A M   O U T P U T
V=NON_VERBOSE################################# Possibilities: 'VERBOSE' / 'NON_VERBOSE'
FILE=FILE_REPRESENTATION###################### Possibilities: 'FILE_REPRESENTATION' / 'NON_FILE_REPRESENTATION'
CPG=CPGPLOT_REPRESENTATION#################### Possibilities: 'CPGPLOT_REPRESENTATION' / 'NON_CPGPLOT_REPRESENTATION'
STATIONARITY=STATIONARY_POINT_REPRESENTATION## Possibilitites:'STATIONARY_POINT_REPRESENTATION' / 
##############################################                'NON_STATIONARY_POINT_REPRESENTATION'
##############################################################################################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   P A R A M E T E R   S P A C E ###############
TYPE_of_BOUNDARY=###
TYPE_of_INDEX=_5p#   '_X2-Y_5p' /
##############################################################################################################
RANDOM_SEED=CURRENT_TIME_RDN_SEED# Alternative: 'FILE_RND_SIDE'
##############################################################################################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   M O D E L   D E F I N I T I O N #############
MODEL=DIFFUSION_S_RESOURCES#  ######################################
ODE_MODEL=${MODEL}####################
################################################### (see rainTemperature.c) ##################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   N U M E R I C A L   I N T E G R A T I O N  ##
ODE_DIR_GENERAL=${PROJECT_DIR}Definition_Numerical_Integration/#########
ODE_DIR=${ODE_DIR_GENERAL}${ODE_MODEL}/#########
ODE_DEF=${ODE_DIR_GENERAL}ODE_Definitions/###########################
ODE_FIXED_POINTS=${PROJECT_DIR}Definition_Fixed_Points/${ODE_MODEL}/###########
TYPE_NUMERICAL_INTEGRATION=GSL_RKF45#############################
######################################## Possibilities:
######################################## 'GSL_BSIMP' uses the Jacobian matrix which has to be 
########################################  provided in JAC_sys${SYS}.c.
######################################## 'GSL_RKF45' is a general purpose integrator.
######################################## 'GSL_BSIMP' and GSL_RKF45 uses adaptive steps
######################################## 'GSL_RK4' does not uses adaptive steps
##############################################################################################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   J A C O B I A N   T Y P E   #################
## To pass this 'JACOBIAN' variable to the compiling code is necessary whether you need to calculate stationary
## points or perform the numerical integration of the system
JACOBIAN=ODE# Options: 'ODE'
EXTRA_INCLUDE_PATH+=-I${ODE_DEF}Include_${JACOBIAN}/#####
###################################o###########################################################################
################################################### (see rainTemperature.c) ###################################
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   S T O C H A S T I C   R E A L I Z A T I O N S#
STO_REA=STO_REALIZATIONS################ Possibilities: 'NON_STO_REALIZATIONS' or 'STO_REALIZATIONS'
STO_DIR_GENERAL=${PROJECT_DIR}Definition_Stochastic_Realizations/#########
STO_DIR=${STO_DIR_GENERAL}${ODE_MODEL}/#########
## D E F I N I T I O N   N A M E   E X E C U T B L E   C O D E ###############################################
EXE=${ODE_MODEL}
##############################################################################################################
EXTERN=#Possibilities: void (to be locally executed) / -static (to be non-locally executed) #
######################################## Possibilities: GA (locally executed) / GA_EXTERN (non-locally executed)
LDFLAGS+=-lm -lgsl -lgslcblas -lDa_GSL_UTL -lda_IO -lda_Optimization -lda_Error -lNumRec -lNRUTL  -lda_cpgplot_XY_GRID_FILE -lda_cpgplot_BASIC ${EXTRA_LIBRARY_PATH}
CXXFLAGS+=${EXTRA_INCLUDE_PATH} -D${V} -D${CPG} -D${FILE}  -D${TYPE_NUMERICAL_INTEGRATION} -D${STO_REA} -D${RANDOM_SEED} -D${MODEL} -D${JACOBIAN} -D${NR} -D${STATIONARITY} -D${STO} -D${BT} -D${RE} -DHAVE_INLINE ${EXTERN} -g #-Wall -D${ERROR}
##############################################################################################################
# either -g or some optimization opiton, for instance, -O3
# This makefile uses full power of conditional compilation. There is a bunch of variables that control this and
# are passed as flags -D[NAME_of_the_VARIABLE] to the compiler during compilation
##############################################################################################################
# Mac users, please do: 
# LIBS=-L/opt/X11/lib -lX11.6
# PGPLOT_LIB=-L/opt/local/lib -lpgplot 
# CPGPLOT_LIB=-L/opt/local/lib -lcpgplot -lpgplot -L/opt/local/lib -lpng
# Linux users please do this insteadd:
LIBS=-L/usr/lib -lX11
PGPLOT_LIB=-L/usr/lib -lpgplot -L/usr/lib -lpng
CPGPLOT_LIB=-L/usr/lib -lcpgplot -lpgplot -L/usr/lib -lpng
##############################################################################################################
OBJECTS = main.o report.o argumentControl.o MODEL.o MODEL_STO.o Community.o Community_Scatter_Plot_Representation.o           \
	Parameter_Model.o Time_Control.o Parameter_Table.o Setting_up_Constant_Metacommunity_Matrix.o                             \
	assignCodes_to_Model_Parameters.o assignLabel_to_Model_Parameters.o                                                       \
	assignLabel_to_Model_Variables.o assignLabel_to_Output_Variables.o                                                        \
	assignStructValue_to_VectorEntry.o assignVectorEntry_to_Structure.o                                                       \
	Vector_Entry_into_Model_Variable_Initial_Condition.o Model_Variable_Initial_Condition_into_Vector_Entry.o                 \
	assert_right_model_definition.o                                                                                           \
	assignSymbol_to_Model_Parameters.o assignSymbol_to_Model_Variables.o assignSymbol_to_Output_Variables.o                   \
	assignLatex_Label_to_Model_Parameters.o  Latex_Write_Parameter_Table.o                                                    \
	fprintf_Model_Parameters.o fprintf_Output_Variables.o                                                                     \
	Parameter_Space.o Parameter_Space_Boundary_for_Trends${TYPE_of_BOUNDARY}.o  Initial_Condition_Boundary_Values.o           \
	deterministic_time_dynamics.o numerical_Integration_Driver.o                                                              \
	Fixed_Points_All.o Stationary_Solution_dydt_Check.o Stationary_Solution_Feasibility_Control.o Stability_2D_Systems.o      \
	Function_Type_of_Stability.o                                                                                              \
	Initial_Conditions_Numerical_Integration.o Initial_Condition_from_Parameter_Table.o                                       \
	Random_Initial_Condition.o ODE_sys_${ODE_MODEL}.o JAC_sys_MODEL.o JACOBIAN_Matrix_${ODE_MODEL}.o                          \
	Stochastic_Time_Dynamics.o Initial_Conditions_Stochastic_Dynamics.o Patch_System_Initialization.o                         \
	Choose_Village.o Advance_Current_Time.o Temporal_Dynamics.o Temporal_Dynamics_Update.o Execute_One_Step.o                 \
	Event_Delta_Matrix_Initialization.o Event_Adjacence_List_Initialization.o Temporal_Dynamics_Aux_Functions.o               \
	Model_Variables_Code.o definition_OutPut_Variables.o definition_OutPut_Variables_Functions.o In_Out_Migration_Functions.o \
	stat.o                                                                                                                    \
	CPGPLOT___X_Y___Parameter_Table.o                                                                                         \
	CPGPLOT___GRID___Parameter_Table.o                                                                                        \
	Time_Dependence_Apply.o Time_Dependence_Resolve.o Time_Dependence_Control.o Trend_Control.o                               
# # index.Definition${TYPE_of_INDEX}.o Error_Function.o 
### Example:

RUN_ARGS =-y0 1 -y2 1 -HM 100 -HX 10 -HY 10 -Hu 0.1 -n 1 -v0 0 -G0 1 -G1 1 -tn 50 -t0 0.0 -t1 50.0 -t4 0 -tR 10 -xn 0 -xN 1000 \
		   -HN 1000 -G2 1 -G3 0.0 -G4 50.0 -G5 1 -G6 0.0 -G7 6000.0

CC=gcc##
FF=gfortran##

all: ${EXE} COMPILATION_OPTIONS

${EXE}:  	$(OBJECTS)
		$(FF) -o ${EXE} $(OBJECTS) $(LDFLAGS) $(CXXFLAGS) $(CPGPLOT_LIB) $(LIBS)

# Main program functions:
main.o: main.c ${HEAD.H}
	$(CC) -c main.c $(CXXFLAGS)

report.o: report.c ${HEAD.H}
	${CC} -c report.c $(CXXFLAGS)

argumentControl.o: argumentControl.c ${HEAD.H}
	${CC} -c argumentControl.c $(CXXFLAGS)

MODEL_STO.o: MODEL_STO.c ${HEAD.H}
	$(CC) -c MODEL_STO.c $(CXXFLAGS)

MODEL.o: MODEL.c ${HEAD.H}
	$(CC) -c MODEL.c $(CXXFLAGS)

Community.o: ${PROJECT_DIR}Community.c ${HEAD.H}
	$(CC) -c ${PROJECT_DIR}Community.c $(CXXFLAGS)

Community_Scatter_Plot_Representation.o: ${PROJECT_DIR}Community_Scatter_Plot_Representation.c ${HEAD.H}
	$(CC) -c ${PROJECT_DIR}Community_Scatter_Plot_Representation.c $(CXXFLAGS)

# Funcions controling main data structures (Parameter_Model, Time_Control,
# Parameter_Table):
assert_right_model_definition.o: ${PROJECT_DIR}assert_right_model_definition.c ${HEAD.H}
	$(CC) -c ${PROJECT_DIR}assert_right_model_definition.c $(CXXFLAGS)

Parameter_Model.o: ${PROJECT_DIR}Parameter_Model.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Parameter_Model.c $(CXXFLAGS)

Time_Control.o: ${PROJECT_DIR}Time_Control.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Time_Control.c $(CXXFLAGS)

Trend_Control.o: ${TIME_TREND_DIR}Trend_Control.c ${HEAD.H}
	${CC} -c ${TIME_TREND_DIR}Trend_Control.c $(CXXFLAGS)

Time_Dependence_Apply.o: ${TIME_TREND_DIR}Time_Dependence_Apply.c ${HEAD.H}
	${CC} -c ${TIME_TREND_DIR}Time_Dependence_Apply.c $(CXXFLAGS)

Time_Dependence_Resolve.o: ${TIME_TREND_DIR}Time_Dependence_Resolve.c ${HEAD.H}
	${CC} -c ${TIME_TREND_DIR}Time_Dependence_Resolve.c $(CXXFLAGS)

Time_Dependence_Control.o: ${TIME_TREND_DIR}Time_Dependence_Control.c ${HEAD.H}
	${CC} -c ${TIME_TREND_DIR}Time_Dependence_Control.c $(CXXFLAGS)

Setting_up_Constant_Metacommunity_Matrix.o: ${PROJECT_DIR}Setting_up_Constant_Metacommunity_Matrix.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Setting_up_Constant_Metacommunity_Matrix.c $(CXXFLAGS)	

Parameter_Table.o: ${PROJECT_DIR}Parameter_Table.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Parameter_Table.c $(CXXFLAGS)

definition_OutPut_Variables.o: ${PROJECT_DIR}definition_OutPut_Variables.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}definition_OutPut_Variables.c ${CXXFLAGS}

definition_OutPut_Variables_Functions.o: ${PROJECT_DIR}definition_OutPut_Variables_Functions.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}definition_OutPut_Variables_Functions.c ${CXXFLAGS}

assignLabel_to_Model_Parameters.o: ${PROJECT_DIR}assignLabel_to_Model_Parameters.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}assignLabel_to_Model_Parameters.c $(CXXFLAGS)

assignCodes_to_Model_Parameters.o:	${PROJECT_DIR}/assignCodes_to_Model_Parameters.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}$/assignCodes_to_Model_Parameters.c $(CXXFLAGS)

assignSymbol_to_Model_Parameters.o:	${PROJECT_DIR}/assignSymbol_to_Model_Parameters.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}$/assignSymbol_to_Model_Parameters.c $(CXXFLAGS)

assignLabel_to_Model_Variables.o: ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/assignLabel_to_Model_Variables.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/assignLabel_to_Model_Variables.c $(CXXFLAGS)

assignSymbol_to_Model_Variables.o: ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/assignSymbol_to_Model_Variables.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/assignSymbol_to_Model_Variables.c $(CXXFLAGS)

assignLabel_to_Output_Variables.o: ${PROJECT_DIR}assignLabel_to_Output_Variables.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}assignLabel_to_Output_Variables.c $(CXXFLAGS)

assignSymbol_to_Output_Variables.o: ${PROJECT_DIR}assignSymbol_to_Output_Variables.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}assignSymbol_to_Output_Variables.c $(CXXFLAGS)

assignVectorEntry_to_Structure.o: ${PROJECT_DIR}assignVectorEntry_to_Structure.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}assignVectorEntry_to_Structure.c $(CXXFLAGS)

assignStructValue_to_VectorEntry.o: ${PROJECT_DIR}assignStructValue_to_VectorEntry.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}assignStructValue_to_VectorEntry.c $(CXXFLAGS)

assignLatex_Label_to_Model_Parameters.o: ${PROJECT_DIR}assignLatex_Label_to_Model_Parameters.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}assignLatex_Label_to_Model_Parameters.c  ${CXXFLAGS}

Latex_Write_Parameter_Table.o:  ${PROJECT_DIR}/Latex_Write_Parameter_Table.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Latex_Write_Parameter_Table.c ${CXXFLAGS}

Model_Variable_Initial_Condition_into_Vector_Entry.o: ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/Model_Variable_Initial_Condition_into_Vector_Entry.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/Model_Variable_Initial_Condition_into_Vector_Entry.c $(CXXFLAGS)

Vector_Entry_into_Model_Variable_Initial_Condition.o: ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/Vector_Entry_into_Model_Variable_Initial_Condition.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Assign_Model_Variables/${ODE_MODEL}/Vector_Entry_into_Model_Variable_Initial_Condition.c $(CXXFLAGS)

fprintf_Model_Parameters.o: ${PROJECT_DIR}fprintf_Model_Parameters.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}fprintf_Model_Parameters.c $(CXXFLAGS)

fprintf_Output_Variables.o: ${PROJECT_DIR}fprintf_Output_Variables.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}fprintf_Output_Variables.c $(CXXFLAGS)


# Functions to control the definition of a subspace from the whole parameter space
Parameter_Space.o: ${PROJECT_DIR}Definition_Parameter_Space/Parameter_Space.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Definition_Parameter_Space/Parameter_Space.c $(CXXFLAGS)	

Parameter_Space_Boundary_for_Trends${TYPE_of_BOUNDARY}.o:  ${PROJECT_DIR}Definition_Parameter_Space/Parameter_Space_Boundary_for_Trends${TYPE_of_BOUNDARY}.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Definition_Parameter_Space/Parameter_Space_Boundary_for_Trends${TYPE_of_BOUNDARY}.c $(CXXFLAGS)

Initial_Condition_Boundary_Values.o:  ${PROJECT_DIR}Definition_Numerical_Integration/${ODE_MODEL}/Initial_Condition_Boundary_Values.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Definition_Numerical_Integration/${ODE_MODEL}/Initial_Condition_Boundary_Values.c $(CXXFLAGS)

# index.Definition${TYPE_of_INDEX}.o: ${PROJECT_DIR}Definition_Parameter_Space/index.Definition${TYPE_of_INDEX}.c ${HEAD.H}
#	${CC} -c ${PROJECT_DIR}Definition_Parameter_Space/index.Definition${TYPE_of_INDEX}.c $(CXXFLAGS)

# Functions to control graphical plots
CPGPLOT___X_Y___Parameter_Table.o: ${CPGPLOT}CPGPLOT_Parameter_Table/CPGPLOT___X_Y___Parameter_Table.c ${CPGPLOT}CPGPLOT_Parameter_Table/CPGPLOT___X_Y___Parameter_Table.h
	${CC} -c ${CPGPLOT}CPGPLOT_Parameter_Table/CPGPLOT___X_Y___Parameter_Table.c $(CXXFLAGS)

CPGPLOT___GRID___Parameter_Table.o: ${CPGPLOT}CPGPLOT_Parameter_Table/CPGPLOT___GRID___Parameter_Table.c ${CPGPLOT}CPGPLOT_Parameter_Table/CPGPLOT___X_Y___Parameter_Table.h
	${CC} -c ${CPGPLOT}CPGPLOT_Parameter_Table/CPGPLOT___GRID___Parameter_Table.c $(CXXFLAGS)

# Functions controling numerical integration of the ODE system
deterministic_time_dynamics.o: ${ODE_DIR_GENERAL}deterministic_time_dynamics.c ${HEAD.H}
	$(CC) -c ${ODE_DIR_GENERAL}deterministic_time_dynamics.c $(CXXFLAGS)

numerical_Integration_Driver.o: ${ODE_DIR_GENERAL}numerical_Integration_Driver.c ${HEAD.H}
	$(CC) -c ${ODE_DIR_GENERAL}numerical_Integration_Driver.c $(CXXFLAGS)

Initial_Conditions_Numerical_Integration.o: ${ODE_DIR_GENERAL}Initial_Conditions_Numerical_Integration.c ${HEAD.H}
	${CC} -c ${ODE_DIR_GENERAL}Initial_Conditions_Numerical_Integration.c $(CXXFLAGS)

Random_Initial_Condition.o: ${ODE_DIR_GENERAL}Random_Initial_Condition.c ${HEAD.H}
	${CC} -c ${ODE_DIR_GENERAL}Random_Initial_Condition.c $(CXXFLAGS)

Initial_Condition_from_Parameter_Table.o: ${ODE_DIR}Initial_Condition_from_Parameter_Table.c ${HEAD.H}
	${CC} -c ${ODE_DIR}Initial_Condition_from_Parameter_Table.c $(CXXFLAGS)

JAC_sys_MODEL.o: ${ODE_DEF}JAC_sys_MODEL.c ${HEAD.H}
	$(CC) -c ${ODE_DEF}JAC_sys_MODEL.c $(CXXFLAGS)

JACOBIAN_Matrix_${ODE_MODEL}.o: ${ODE_DEF}JACOBIAN_Matrix_${ODE_MODEL}.c ${HEAD.H}
	$(CC) -c ${ODE_DEF}JACOBIAN_Matrix_${ODE_MODEL}.c $(CXXFLAGS)

ODE_sys_${ODE_MODEL}.o: ${ODE_DEF}ODE_sys_${ODE_MODEL}.c ${HEAD.H}
	$(CC) -c ${ODE_DEF}ODE_sys_${ODE_MODEL}.c $(CXXFLAGS)

# Functions controling numerical integration of the ODE system
# Stochastic_Time_Dynamics.o Initial_Conditions_Stochastic_Dynamics.o Patch_System_Initialization.o Advance_Current_Time.o Choose_Village.o Temporal_Dynamics.o Execute_One_Step.o Force_of_Infection_Sto.o
Temporal_Dynamics_Update.o: ${STO_DIR}Temporal_Dynamics_Update.c ${HEAD.H}
	${CC} -c ${STO_DIR}Temporal_Dynamics_Update.c $(CXXFLAGS)

Temporal_Dynamics.o: ${STO_DIR}Temporal_Dynamics.c ${HEAD.H}
	${CC} -c ${STO_DIR}Temporal_Dynamics.c $(CXXFLAGS)

Execute_One_Step.o: ${STO_DIR}Execute_One_Step.c ${HEAD.H}
	${CC} -c ${STO_DIR}Execute_One_Step.c $(CXXFLAGS)

Event_Delta_Matrix_Initialization.o: ${STO_DIR}Event_Delta_Matrix_Initialization.c ${HEAD.H}
	${CC} -c ${STO_DIR}Event_Delta_Matrix_Initialization.c $(CXXFLAGS)

Event_Adjacence_List_Initialization.o: ${STO_DIR}Event_Adjacence_List_Initialization.c ${HEAD.H}
	${CC} -c ${STO_DIR}Event_Adjacence_List_Initialization.c $(CXXFLAGS)

Temporal_Dynamics_Aux_Functions.o: ${STO_DIR}Temporal_Dynamics_Aux_Functions.c ${HEAD.H}
	${CC} -c ${STO_DIR}Temporal_Dynamics_Aux_Functions.c $(CXXFLAGS)

Stochastic_Time_Dynamics.o: ${STO_DIR_GENERAL}Stochastic_Time_Dynamics.c ${HEAD.H}
	${CC} -c ${STO_DIR_GENERAL}Stochastic_Time_Dynamics.c $(CXXFLAGS)

Initial_Conditions_Stochastic_Dynamics.o: ${STO_DIR_GENERAL}Initial_Conditions_Stochastic_Dynamics.c ${HEAD.H}
	${CC} -c ${STO_DIR_GENERAL}Initial_Conditions_Stochastic_Dynamics.c $(CXXFLAGS)

Patch_System_Initialization.o: ${STO_DIR_GENERAL}Patch_System_Initialization.c ${HEAD.H}
	${CC} -c ${STO_DIR_GENERAL}Patch_System_Initialization.c $(CXXFLAGS)

Advance_Current_Time.o: ${STO_DIR_GENERAL}Advance_Current_Time.c ${HEAD.H}
	${CC} -c ${STO_DIR_GENERAL}Advance_Current_Time.c $(CXXFLAGS)

Choose_Village.o: ${STO_DIR_GENERAL}Choose_Village.c ${HEAD.H}
	${CC} -c ${STO_DIR_GENERAL}Choose_Village.c $(CXXFLAGS)

# Auxiliary functions to specify model ODE definition
Model_Variables_Code.o: ${PROJECT_DIR}Model_Variables_Code.c ${HEAD.H}
	${CC} -c ${PROJECT_DIR}Model_Variables_Code.c $(CXXFLAGS)

In_Out_Migration_Functions.o: ${PROJECT_DIR}In_Out_Migration_Functions.c ${HEAD.H}
	$(CC) -c ${PROJECT_DIR}In_Out_Migration_Functions.c  $(CXXFLAGS)

# # # Funcions to calculate the fixed points of the ODE system
Fixed_Points_All.o: ${ODE_FIXED_POINTS}Fixed_Points_All.c ${HEAD.H}
	$(CC) -c ${ODE_FIXED_POINTS}Fixed_Points_All.c $(CXXFLAGS)

Stationary_Solution_dydt_Check.o: ${FIXED_POINTS_DIR}Stationary_Solution_dydt_Check.c ${HEAD.H}
	$(CC) -c ${FIXED_POINTS_DIR}Stationary_Solution_dydt_Check.c $(CXXFLAGS)

Stationary_Solution_Feasibility_Control.o: ${FIXED_POINTS_DIR}Stationary_Solution_Feasibility_Control.c ${HEAD.H}
	$(CC) -c ${FIXED_POINTS_DIR}Stationary_Solution_Feasibility_Control.c $(CXXFLAGS)

Stability_2D_Systems.o: ${FIXED_POINTS_DIR}Stability_2D_Systems.c ${HEAD.H}
	${CC} -c ${FIXED_POINTS_DIR}Stability_2D_Systems.c $(CXXFLAGS)

Function_Type_of_Stability.o: ${FIXED_POINTS_DIR}Function_Type_of_Stability.c ${HEAD.H}
	${CC} -c ${FIXED_POINTS_DIR}Function_Type_of_Stability.c $(CXXFLAGS)
			
# fixed_ALL_Points_Intersection.o: ${ODE_FIXED_POINTS}fixed_ALL_Points_Intersection.c ${HEAD.H}
#	$(CC) -c ${ODE_FIXED_POINTS}fixed_ALL_Points_Intersection.c $(CXXFLAGS)

# Steady_State_x_Fraction.o: ${ODE_FIXED_POINTS}Steady_State_x_Fraction.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}Steady_State_x_Fraction.c $(CXXFLAGS)

# Steady_State_y_Fraction.o: ${ODE_FIXED_POINTS}Steady_State_y_Fraction.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}Steady_State_y_Fraction.c $(CXXFLAGS)

# x_y_Fixed_Points_Functions.o: ${ODE_FIXED_POINTS}x_y_Fixed_Points_Functions.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}x_y_Fixed_Points_Functions.c $(CXXFLAGS)

# y_W_diagramm.o: ${ODE_FIXED_POINTS}y_W_diagramm.c ${HEAD.H}
#	$(CC) -c ${ODE_FIXED_POINTS}y_W_diagramm.c $(CXXFLAGS)

# intermediateRoot_x_Fraction.o: ${ODE_FIXED_POINTS}intermediateRoot_x_Fraction.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}intermediateRoot_x_Fraction.c $(CXXFLAGS)

# Stable_Solution.o: ${ODE_FIXED_POINTS}Stable_Solution.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}Stable_Solution.c $(CXXFLAGS)

# Stable_Solution_x.o: ${ODE_FIXED_POINTS}${X_MODEL}/Stable_Solution_x.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}${X_MODEL}/Stable_Solution_x.c $(CXXFLAGS)

# Stable_Solution_y.o: ${ODE_FIXED_POINTS}${Y_MODEL}/Stable_Solution_y.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}${Y_MODEL}/Stable_Solution_y.c $(CXXFLAGS)

# Generic_Root_Finder.o: ${ODE_FIXED_POINTS}Generic_Root_Finder.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}Generic_Root_Finder.c $(CXXFLAGS)

# Feasibility_Condition.o: ${ODE_FIXED_POINTS}Feasibility_Condition.c ${HEAD.H}
# 	$(CC) -c ${ODE_FIXED_POINTS}Feasibility_Condition.c $(CXXFLAGS)

stat.o:	${PROJECT_DIR}stat.c ${PROJECT_DIR}stat.h
	${CC} -c ${PROJECT_DIR}stat.c $(CXXFLAGS)
#

clean:
	rm -f *.o *~ core
#
COMPILATION_OPTIONS:
	echo "Main model definition compilation options:"                        > ${EXE}_COMPILATION_OPTIONS.txt
	echo "ODE algorithm: ${TYPE_NUMERICAL_INTEGRATION}"                      >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "Graphical representation: ${CPG}"                                  >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "Stochastic optimization: "
	echo "Binary tree optimization: ${BT}"                                   >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "Random number reuse: ${RE}"                                        >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "Stochastic Optimization: ${STO}"                                   >> ${EXE}_COMPILATION_OPTIONS.txt	 	
	echo "Fixed Points Calculation: ${STATIONARITY}"                         >> ${EXE}_COMPILATION_OPTIONS.txt	
	echo "Stochastic Realizations: ${STO_REA}"                               >> ${EXE}_COMPILATION_OPTIONS.txt	
	echo ""           			              	                             >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "---------------------------------------" >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "All compilation options as passed to the compiler through CXXFLAGS:" >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "${CXXFLAGS}"                             >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "---------------------------------------" >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "All linking options as passed to the linker through LDFLAGS, CPGPLOT_LIG and LIBS:" >> ${EXE}_COMPILATION_OPTIONS.txt
	echo "$(LDFLAGS) $(CPGPLOT_LIB) $(LIBS)"                                                  >> ${EXE}_COMPILATION_OPTIONS.txt

run:
	${WKD}${EXE} ${RUN_ARGS}
