# =======================
#     LIBRARY make
# =======================

# Gnu C make cpgplot_GRID UTILS (libda_cpgplot_GRID.a)
NR=NUMERICAL_RECIPES
CPGPLOT=${HOME}/CPGPLOT/
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
NR_DIR               = ${PROJECT_DIR}Library/C_NumRecip/
LIB_DIR              = ${PROJECT_DIR}Library/Optimization_Library/
RANDOM_DIR           = ${PROJECT_DIR}Library/GSL_RANDOM/
IO_DIR               = ${PROJECT_DIR}Library/IO_Library/
ERROR_DIR            = ${PROJECT_DIR}Definition_Error_Model/
# Linux users, do:
INCLUDE_R_DIR        =/usr/share/R/include/
# Mac urers, please do:
# INCLUDE_DIR        =/opt/local/include/
# Linux users, do:
# EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${PROJECT_DIR}Include/ -I${RANDOM_DIR} -I${IO_DIR} -I${INCLUDE_R_DIR} -I${ERROR_DIR} -I${NR_DIR}
# Mac user, plese do:
# EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${PROJECT_DIR}Include/ -I${RANDOM_DIR} -I${IO_DIR} -I${INCLUDE_DIR} -I${ERROR_DIR}
EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${PROJECT_DIR}Include/ -I${RANDOM_DIR} -I${IO_DIR} -I${INCLUDE_R_DIR} -I${ERROR_DIR} -I${NR_DIR}
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}###################
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library/###########
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library_basic_UTL/#

FILE=FILE_REPRESENTATION###################### Possibilities: 'FILE_REPRESENTATION' / 'NON_FILE_REPRESENTATION'
CPG=CPGPLOT_REPRESENTATION############## Possibilities: 'CPGPLOT_REPRESENTATION' / 'NON_CPGPLOT_REPRESENTATION'
## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   P R O G R A M   O U T P U T
V=NON_VERBOSE###################################### Possibilities: 'VERBOSE' / 'NON_VERBOSE'

CFLAGS=-g -Wall -DHAVE_INLINE -D${CPG} -D${V} ${EXTRA_INCLUDE_PATH} -D${NR}
###CFLAGS= -Wall -O21

COMPILER=gcc -fPIC -DSHARED_OBJECT

OBJECTS=Function_to_Minimize.o generic_Function_Parameter_2Dim_Scan.o generic_Function_Parameter_2Dim_Scan_Improved.o        \
		GSL_Function_to_Minimize.o GSL_Function_to_Minimize_Error_Model.o GSL_Function_to_Minimize_Binomial_Free_Consumers.o \
		GSL_Function_to_Minimize_Multinomial_Free_Consumers.o GSL_Function_to_Minimize_Beddington_DeAngelis.o                \
		GSL_Function_to_Minimize_Beddington_DeAngelis_Marginal_0.o GSL_neglog_Error_Probability_Model.o                      \
		GSL_neglog_Error_Probability_Model_Gaussian.o GSL_neglog_Error_Probability_Model_Gamma.o GSL_Minimization_Simplex.o  \
		GSL_CPGPLOT_Minimization_Simplex.o GSL_Minimization_Driver.o Checking_for_Parameter_Boundaries.o                     \
		Inspecting_Likelihood_of_Final_Solution.o Inspecting_Solution_Driver.o                                               \
		GSL_Eigenvalue_Calculation.o NR_Eigenvalue_Calculation.o Confidence_Intervals_from_Likelihood_Profile.o treenode.o 

.SUFIXES:       .o .c

.c.o:
	${COMPILER} $(CFLAGS) -c $<

da_Optimization: ${OBJECTS}
	ar -q libda_Optimization.a ${OBJECTS}

da_Optimization_SO: ${OBJECTS}
	gcc -shared -o libda_Optimization.so ${OBJECTS}

all: da_Optimization 

clean:
	rm -f *.o *~
