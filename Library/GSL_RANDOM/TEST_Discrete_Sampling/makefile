### WKD =/home/dalonso/PROJECT_STOCHASTIC_DIFFUSION/Library/GSL_RANDOM/TEST_Discrete_Sampling_so/#
### $PWD

# Adding definitions for the path where include files are searched for  
EXTRA_INCLUDE_PATH+=-I${PWD}####################

# Adding definitions for the path where libraries (lib*.a) are searched for  
EXTRA_LIBRARY_PATH+=-L${PWD}####################

# Adding definitions for the path where the system looks for shared objects at running time
# :~$ echo $LD_LIBRARY_PATH
# :~$ LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}
# :~$ export LD_LIBRARY_PATH
# After generating the executable, check with 
# :~$ ldd program
# You can do this in a perennial way by editing the .bashrc file, and adding these two lines
# at the end of it:
# LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}
# export LD_LIBRARY_PATH

OBJ = main.o # Discrete_Sampling.o 
LIBRARY_OBJECTS = Discrete_Sampling.o
              
LDFLAGS+=-lm ${EXTRA_LIBRARY_PATH} -lDA_Discrete_Sampling 
CXXFLAGS+=${EXTRA_INCLUDE_PATH} -g 

CC  = gcc

program:	${OBJ} # shared_object
		${CC} -o $@ ${OBJ} ${LDFLAGS} ${CXXFLAGS}

library: ${LIBRARY_OBJECTS} Discrete_Sampling.h
	ar -q libDA_Discrete_Sampling.a ${LIBRARY_OBJECTS}

shared_object: ${LIBRARY_OBJECTS}
	gcc -shared -o libDA_Discrete_Sampling.so ${LIBRARY_OBJECTS}

clean:	
	rm *.o

cleanall: 
	rm *.o program libDA_Discrete_Sampling.a libDA_Discrete_Sampling.so