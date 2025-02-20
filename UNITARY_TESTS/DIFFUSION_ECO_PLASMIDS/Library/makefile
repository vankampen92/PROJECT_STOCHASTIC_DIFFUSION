# =======================
#     LIBRARY make
# =======================

# Gnu C make cpgplot_GRID UTILS (libda_cpgplot_GRID.a)
CPGPLOT=${HOME}/CPGPLOT/
PROJECT_DIR=${HOME}/PROJECT_STOCHASTIC_DIFFUSION/
LIB_DIR              = ${PROJECT_DIR}UNITARY_TESTS/DIFFUSION_ECO_PLASMIDS/Library/
INCLUDE_DIR          = ${PROJECT_DIR}UNITARY_TESTS/DIFFUSION_ECO_PLASMIDS/Include/
# Linux users, do:
INCLUDE_R_DIR        =/usr/share/R/include/
# Mac urers, please do:
# INCLUDE_DIR       +=/opt/local/include/

EXTRA_INCLUDE_PATH   = -I${LIB_DIR} -I${PROJECT_DIR} -I${INCLUDE_DIR} 

# EXTRA_INCLUDE_PATH+=-I${CPGPLOT}###################
# EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library/###########
# EXTRA_INCLUDE_PATH+=-I${CPGPLOT}library_basic_UTL/#

# FILE=FILE_REPRESENTATION###################### Possibilities: 'FILE_REPRESENTATION' / 'NON_FILE_REPRESENTATION'
# CPG=CPGPLOT_REPRESENTATION############## Possibilities: 'CPGPLOT_REPRESENTATION' / 'NON_CPGPLOT_REPRESENTATION'

## C O M P I L A T I O N   O P T I O N S   C O N T R O L I N G   P R O G R A M   O U T P U T
V=NON_VERBOSE###################################### Possibilities: 'VERBOSE' / 'NON_VERBOSE'

CFLAGS=-g -Wall -DHAVE_INLINE -D${V} ${EXTRA_INCLUDE_PATH}
###CFLAGS= -Wall -O21

COMPILER=gcc -fPIC -DSHARED_OBJECT

OBJECTS=Calculate_Strain_and_Profile.o \
		Calculate_Strain_ID_from_Profile_and_Strain_Type.o \
		Compatibility_Plasmic_Sets_Assert.o  \
		Compatibility_Plasmids_Sets_Count.o  \
		Compatibility_Profiles_Assert.o  \
		Create_Binary_Combination.o  \
		Determining_actual_No_of_RESOURCES.o  \
		Effective_Donor_Recipient_Reaction_Assert.o  \
		GSL_Init_Random_Seed.o  \
		Infection_Condition_Assert.o  \
		Potential_Donor_Recipient_Pair_Assert.o  \
		Print_Infection_Profile.o  \
		Printing_Putative_Recipient_Lists.o  \
		Printing_Strains_Profiles_and_Lists.o  \
		Printing_Strains_Profiles.o  \
		Print_Strain_Profile.o  \
		Profile_Selfconsistency_Assert.o  \
		Recipient_Donor_Transconjugant_Rate.o  \
		Setting_Adjacency_Lists_from_Interaction_Matrices.o  \
		Setting_Interacction_Matrices.o  \
		Setting_Plasmid_Characteristic_Parameters.o  \
		Setting_Putative_Recipient_Lists_of_Potential_Transconjugants.o  \
		Setting_Reactive_Recipient_Donor_Pairs_and_Rates.o  \
		Setting_Strain_Characteristic_Parameters.o  \
		show_DoubleMatrix.o  \
		SumUP_Profile.o  \
		Tranconjugation_Gain_and_Loss_Total_Rates.o

.SUFIXES:       .o .c

.o.c:
	${COMPILER} $(CFLAGS) -c $<

eco_Plasmids: ${OBJECTS}
	ar -q libeco_Plasmids.a ${OBJECTS}

eco_Plasmids_SO: ${OBJECTS}
	gcc -shared -o libeco_Plasmids.so ${OBJECTS}

all: eco_Plasmids 

clean:
	rm -f *.o *~
