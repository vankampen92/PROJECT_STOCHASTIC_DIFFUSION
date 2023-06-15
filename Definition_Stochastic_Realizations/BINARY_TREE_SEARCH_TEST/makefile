OBJ = main_Tree_Discrete_Distribution.o treenode.o 
LIB = -lm
CC  = gcc
CXXFLAGS+= -g

program:	${OBJ} ${LIB}
		${CC} -o $@ ${OBJ} ${LIB} ${CXXFLAGS}

clean:	
		rm *.o

cleanall: 
		rm *.o program