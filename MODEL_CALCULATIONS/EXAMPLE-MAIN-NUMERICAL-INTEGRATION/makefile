# Relevant directories are saved in ENVIRONMENT VARIABLES
CPGPLOT=${HOME}/CPGPLOT/
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
GSL_DA=${PROJECT_DIR}Library/GSL_RANDOM/##########
WKD=${PROJECT_DIR}MODEL_CALCULATIONS/EXAMPLE-MAIN-NUMERICAL-INTEGRATION/

# Adding definitions for the path where include files are searched for  
EXTRA_INCLUDE_PATH+=-I${GSL_DA}####################
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}###################
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library/###########
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library_basic_UTL/#
EXTRA_INCLUDE_PATH+=-I${WKD}####################### 

# Adding definitions for the path where libraries (lib*.a) are searched for  
EXTRA_LIBRARY_PATH+=-L${GSL_DA}####################
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}###################
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}library/###########
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}library_basic_UTL/#

# Definition of ENVIRONMENT VARIABLES to be passed at the compilation step
# which will control conditional compilation
CPG = CPGPLOT_REPRESENTATION

LDFLAGS+=-lm -lgsl -lgslcblas -lDa_GSL_UTL ${EXTRA_LIBRARY_PATH} -lda_cpgplot_XY_GRID_FILE -lda_cpgplot_BASIC
CXXFLAGS+=${EXTRA_INCLUDE_PATH} -g -D${CPG} #-Wall -D${ERROR}

LIBS=-L/usr/lib -lX11
CPGPLOT_LIB=-L/usr/lib -lcpgplot -lpgplot -L/usr/lib -lpng

# GSL_RNG_SEED=104871
EXE     =  cotton
OBJECTS =  main.o argumentControl.o numerical_Integration_Driver.o JAC_sys_MODEL.o ODE_sys_MODEL.o #Ant2.o 

CC=gcc##

# Linking main with libraries 
random:  	$(OBJECTS) random-library
		$(CC) -o ${EXE} $(OBJECTS) $(LDFLAGS) $(CXXFLAGS) $(CPGPLOT_LIB) $(LIBS)

# Compilation of the library
random-library:	
		cd ${GSL_DA}; make clean; make 

# Compiling main program function:
main.o: main.c 
	$(CC) -c main.c $(CXXFLAGS)

argumentControl.o: argumentControl.c 
	$(CC) -c argumentControl.c $(CXXFLAGS)

numerical_Integration_Driver.o: numerical_Integration_Driver.c 
	$(CC) -c numerical_Integration_Driver.c $(CXXFLAGS)

JAC_sys_MODEL.o: JAC_sys_MODEL.c 
	$(CC) -c JAC_sys_MODEL.c $(CXXFLAGS)

ODE_sys_MODEL.o: ODE_sys_MODEL.c 
	$(CC) -c ODE_sys_MODEL.c $(CXXFLAGS)

# Deleting object files
clean:
	rm -f *.o *~ ${EXE}
