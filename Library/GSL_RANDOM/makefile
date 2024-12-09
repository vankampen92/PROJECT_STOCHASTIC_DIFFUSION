# ==========================================
# CUSTUMIZED GSL-based UTILITY LIBRARY make
# ==========================================
#
# ==========================================
# Gnu C make to build GSL-based UTILITY Library 
# The generated library is libDa_GSL_UTL.a 
# ==========================================
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
GSL_Stat_UTL=$(PROJECT_DIR)/Library/GSL_RANDOM/
EXTRA_INCLUDE_PATH= -I${GSL_Stat_UTL}
CFLAGS=-g -fPIE -DHAVE_INLINE ${EXTRA_INCLUDE_PATH}
###CFLAGS= -Wall -O21 
COMPILER=gcc

OBJECTS = GSL_Discrete_Sampling.o GSL_Init_Random_Seed.o GSL_stat.o Percentile_Calculation_Simple.o Create_Binary_Combination.o 

.c.o: 
	${COMPILER} $(CFLAGS) -c $<

GSL_UTL: ${OBJECTS} GSL_stat.h
	ar -q libDa_GSL_UTL.a ${OBJECTS}

all: GSL_UTL

clean:  
	rm -f *.o *~ *.a
