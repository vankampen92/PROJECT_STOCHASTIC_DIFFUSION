#include "../../Include/MODEL.h"
/* Global Shared parameters main Program <---> ArgumentControl() */

#include "extern.h"

void ArgumentControl_fprintf(FILE * fp, Parameter_Table * Table,
			     int argc, char **argv)
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
	      #include <include.Trend_Control.argumentControl_fprintf.c>
	      #include <include.Output_Variables.argumentControl_fprintf.c>
	      #include <include.Parameter_Model.argumentControl_fprintf.c>
	      #include <include.Parameter_Space.argumentControl_fprintf.c>
	      #include <include.Time_Control.argumentControl_fprintf.c>
	      #include <include.Time_Dependence_Control.argumentControl_fprintf.c>

              #include <include.Initial_Conditions.argumentControl_fprintf.c>

	      #include <include.Error_Control.argumentControl_fprintf.c>

  	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    
	      printf("\n");
	      printf(" As an example,\n");
	      printf(" ~$ %s -n 9 -v0 0 -v1 1 -v2 2 -v3 6 -v4 7 -v5 8 -v6 9 -v7 10 -v8 12 -G0 3 -G1 3 -tn 5000 -t0 0.0 -t1 50.0 -t4 0 -xn 0 -xN 1000000 -H2 25000 -H12 25000\n\n", argv[0]);

	      exit(0);
	    }
	}
      else
	{
	  printf("**No command line arguments are written or saved!!!\n");
	  printf("**Invalid command line flag marker (the label for a input argument always starts with a dash) >%c<\n",argv[argcount][0]);
	  
	}
    }
}
