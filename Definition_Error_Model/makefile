# =======================
#     LIBRARY make
# =======================

# Gnu C make:
# Gnu C make cpgplot_GRID UTILS (libda_cpgplot_GRID.a)
CPGPLOT=${HOME}/CPGPLOT/
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
LIB_DIR              = ${PROJECT_DIR}Library/Optimization_Library/
RANDOM_DIR           = ${PROJECT_DIR}Library/GSL_RANDOM/
IO_DIR               = ${PROJECT_DIR}Library/IO_Library/
ERROR_DIR            = ${PROJECT_DIR}Definition_Error_Model/
# Linux users do:
INCLUDE_R_DIR        =/usr/share/R/include/
# Mac urers, please do:
# INCLUDE_DIR        =/opt/local/include/
EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${PROJECT_DIR}Include/ -I${RANDOM_DIR} -I${IO_DIR} -I${INCLUDE_R_DIR} -I${ERROR_DIR}
# Mac users, please do:
# EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${PROJECT_DIR}Include/ -I${RANDOM_DIR} -I${IO_DIR} -I${INCLUDE_DIR} -I${ERROR_DIR}
# Linux users, do this instead: 
EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${PROJECT_DIR}Include/ -I${RANDOM_DIR} -I${IO_DIR} -I${INCLUDE_R_DIR} -I${ERROR_DIR}
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}###################
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library/###########
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library_basic_UTL/#

FILE=FILE_REPRESENTATION###################### Possibilities: 'FILE_REPRESENTATION' / 'NON_FILE_REPRESENTATION'
CPG=CPGPLOT_REPRESENTATION############## Possibilities: 'CPGPLOT_REPRESENTATION' / 'NON_CPGPLOT_REPRESENTATION'

CFLAGS=-g -Wall -DHAVE_INLINE -D${CPG} ${EXTRA_INCLUDE_PATH}
###CFLAGS= -Wall -O21

COMPILER=gcc -fPIC -DSHARED_OBJECT

OBJECTS=Error_Function.o Error_Model_Structure_Value_into_Vector_Entry.o Vector_Entry_into_Error_Model_Structure.o Random_Error_Model_within_Boundaries.o

.SUFIXES:       .o .c

.c.o:
	${COMPILER} $(CFLAGS) -c $<

da_Error: ${OBJECTS}
	ar -q libda_Error.a ${OBJECTS}

da_Error_SO: ${OBJECTS}
	gcc -shared -o libda_Error.so ${OBJECTS}

clean:
	rm -f *.o *~
