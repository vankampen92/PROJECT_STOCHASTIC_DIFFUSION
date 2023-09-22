OBJ = main_Priority_Queu.o treenode.o 
LIB = -lm
CC  = gcc -g
CXXFLAGS+=

program:	${OBJ} ${LIB}
		${CC} -o $@ ${OBJ} ${LIB} ${CXXFLAGS}

clean:	
		rm *.o

cleanall: 
		rm *.o program
