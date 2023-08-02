#include <MODEL.h>
/* Global Shared parameters main Program <---> ArgumentControl() */

#include "extern.h"

extern int * Index_Output_Variables; 

extern int HELP_INPUT_ARGUMENTS;

void ArgumentControl(int argc, char **argv)
{
  int argcount, skip;
#if defined FILE_REPRESENTATION
  int n_Dummy;
#endif

  for(argcount=1; argcount<argc; argcount+=skip)
    {
      if(argv[argcount][0] == '-')
	{
	  skip = 1;
	  switch(argv[argcount][1])
	    {

#if defined CPGPLOT_REPRESENTATION
#include <include.CPG.argumentControl.c>
#if defined FILE_REPRESENTATION
#include <include.FILES_to_READ.argumentControl.c>
#endif
#endif

#include <include.Trend_Control.argumentControl.c>
#include <include.Output_Variables.argumentControl.c>
#include <include.Parameter_Model.argumentControl.c>
#include <include.Parameter_Space.argumentControl.c>
#include <include.Time_Control.argumentControl.c>
#include <include.Time_Dependence_Control.argumentControl.c>
#include <include.Initial_Conditions.argumentControl.c>
#include <include.Error_Control.argumentControl.c>

  	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':

	      P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C( P_ARG );
	      P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( P_ARG, Index_Output_Variables );
	      Parameters_from_Command_Line(stdout, P_ARG);

	      printf("\n");
	      printf(" As an example,\n");
	      printf("   ~$ %s -sP 5 -n 9 -v0 0 -v1 1 -v2 2 -v3 6 -v4 7 -v5 8 -v6 9 -v7 10 -v8 12 -G0 3 -G1 3 -tn 5000 -t0 0.0 -t1 50.0 -t4 0 -xn 0 -xN 1000000 -H2 25000 -H12 25000\n\n", argv[0]);

	      P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( P_ARG );
	      exit(0);
	    }
	}
      else
	{
	  printf("**invalid command line flag >%c<\n",argv[argcount][0]);
	  printf("try -h for help.\n");
	  exit(0);
	}
    }
}
