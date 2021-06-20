#include <MODEL.h>
/* Global Shared parameters main Program <---> ArgumentControl() */

#include "extern.h"

void ArgumentControl_fprintf(FILE * fp, Parameter_Table * Table, int argc, char **argv)
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
#include <include.CPG.argumentControl_fprintf.c>
#if defined FILE_REPRESENTATION
#include <include.FILES_to_READ.argumentControl_fprintf.c>
#endif
#endif
	      
/* #include <include.Trend_Control.argumentControl_fprintf.c> */
/* #include <include.Output_Variables.argumentControl_fprintf.c> */
/* #include <include.Parameter_Model.argumentControl_fprintf.c> */
/* #include <include.Parameter_Space.argumentControl_fprintf.c>  */
/* #include <include.Time_Control.argumentControl_fprintf.c>  */
/* #include <include.Time_Dependence_Control.argumentControl_fprintf.c>  */
#include <include.Initial_Conditions.argumentControl_fprintf.c>
#include <include.Error_Control.argumentControl_fprintf.c>  
	      
	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	      
	      printf("\n");
	      printf(" As an example,\n");
	      printf(" ~$ %s -y0 1 -G0 2 -G1 2 -n 4 -v0 3 -v1 4 -v2 5 -v3 6 -sT 1.0E-04 -sN 300 -sP 2 -I0 21 -H21 1.0 -m0 0.8 -M0 1.2 -A0 0.01 -I1 0 -H0 100.0 -m1 96.0 -M1 120.0 -A1 0.1 -iP 0 -en 0 -eV 10.0  -tn 100 -t0 0.0 -t1 10.0 -t4 1 -xn 0 -xN 0.0 -tE 2.0 -tR 100 -xn 0 -xR 1 -xN 0.0 -DP 1 -DC- -D0 0 -D1 1 -D2 0 -P0 16 -a0 0 -Fn 2 -F0 Pseudo_Data_File.dat -F1 Time_Dependent_Downloading_Rate.dat \n\n", argv[0]);
	      
	      exit(0);
	    }  
	}
      else
	{
	  printf("**No command line arguments will be written or saved!!!\n");
	  printf("**Invalid command line flag marker (the label for a input argument always starts with a dash) >%c<\n", argv[argcount][0]);
	}
    }
}

