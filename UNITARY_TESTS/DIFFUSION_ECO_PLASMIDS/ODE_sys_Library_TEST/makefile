
CPGPLOT=${HOME}/CPGPLOT/
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/UNITARY_TESTS/DIFFUSION_ECO_PLASMIDS/ODE_sys_TEST/

EXTRA_LIBRARY_PATH+=-L${PROJECT_DIR}###############
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}library_basic_UTL/
EXTRA_LIBRARY_PATH+=-L${CPGPLOT}library/

EXTRA_INCLUDE_PATH+=-I${PROJECT_DIR}###############
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library/
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library_basic_UTL/
EXTRA_INCLUDE_PATH+=-I${CPGPLOT}

CPG=CPGPLOT_REPRESENTATION#############

OBJ = main.o argumentControl.o functions.o deterministic_time_dynamics.o Initial_Conditions_Numerical_Integration.o numerical_Integration_Driver.o Parameter_Table.o ODE_sys_DIFFUSION_ECO_PLASMIDS.o    
LIB = -lm -lgsl -lgslcblas ${EXTRA_LIBRARY_PATH} -lda_cpgplot_XY_GRID_FILE -lda_cpgplot_BASIC
PLPLOT_LIB=-L/usr/lib -lcpgplot -lpgplot -lX11 -lpng

CXXFLAGS+= ${EXTRA_INCLUDE_PATH} -D${CPG}
CC  = gcc -g ${CXXFLAGS} 

program:	${OBJ} ${LIB}
			${CC} -o $@ ${OBJ} ${LIB}

plot: ${OBJ}
			gfortran $(CXXFLAGS) -o plot $(OBJ) $(LIB) $(PLPLOT_LIB)

clean:	
		rm *.o

cleanall: 
		rm *.o program
