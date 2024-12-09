OBJ = main.o  
LIB = -lm -lgsl -lgslcblas
CC  = gcc -g
CXXFLAGS+=

program:	${OBJ} ${LIB}
		${CC} -o $@ ${OBJ} ${LIB} ${CXXFLAGS}

clean:	
		rm *.o

cleanall: 
		rm *.o program

