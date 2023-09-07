OBJ = main_Tree_Discrete_Distribution_Bis.o treenode.o 
LIB = -lm
CC  = gcc -g
CXXFLAGS+=

program:	${OBJ} ${LIB}
		${CC} -o $@ ${OBJ} ${LIB} ${CXXFLAGS}

clean:	
		rm *.o

cleanall: 
		rm *.o program
