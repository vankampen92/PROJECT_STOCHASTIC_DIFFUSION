PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
GSL_DA=${PROJECT_DIR}Library/GSL_RANDOM/##########
EXTRA_INCLUDE_PATH+=-I${GSL_DA}####################
EXTRA_LIBRARY_PATH+=-L${GSL_DA}####################

LDFLAGS+=${EXTRA_LIBRARY_PATH} -lm -lgsl -lgslcblas -lDa_GSL_UTL 
CXXFLAGS+=${EXTRA_INCLUDE_PATH} -g 

OBJECTS =  main.o

CC=gcc##

# Linking main with libraries 
random:  	$(OBJECTS) random-library
		$(CC) -o random $(OBJECTS) $(LDFLAGS) $(CXXFLAGS)

# Compilation of the library
random-library:	
		cd ${GSL_DA}; make clean; make 

# Compiling main program function:
main.o: main.c 
	$(CC) -c main.c $(CXXFLAGS)

# Deleting object files
clean:
	rm -f *.o *~ random
