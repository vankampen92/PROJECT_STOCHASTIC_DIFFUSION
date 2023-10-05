OBJ = main_Binary_Tree.o treenode.o 
LIB = -lm
CC  = gcc -g
CXXFLAGS+=

program:	${OBJ} ${LIB}
		${CC} -o $@ ${OBJ} ${LIB} ${CXXFLAGS}

clean:	
		rm *.o

cleanall: 
		rm *.o program
