#include <MODEL.h>
/* Global Shared parameters main Program <---> ArgumentControl() */

#include "extern.h"

void ArgumentControl(int argc, char **argv)
{
  int argcount, skip;
  
  int n_Dummy;  /* FILE_REPRESENTATION dummy variable */ 
  
  for(argcount=1; argcount<argc; argcount+=skip)
    {
      if(argv[argcount][0] == '-')
	{
	  skip = 1;
	  switch(argv[argcount][1])
	    {

#if defined CPGPLOT_REPRESENTATION 
#include <include.CPG.argumentControl.c>
#include <include.FILES_to_READ.argumentControl.c>
#endif

#include <include.Parameter_Model.argumentControl.c>

  	    default:
	      printf("**invalid command line argument >%c< \n",
		     argv[argcount][1]);
	    case 'h':

	      printf("\n");
	      printf(" An example for a command line program call is:\n");
	      printf(" ~$ %s -pK [VALUE] -pBC [VALUE] -pBR [VALUE] -pU [VALUE] -pA [VALUE] -pDC [VALUE] -pDR [VALUE]\n", argv[0]);

	      exit(0);
	    }
	}
      else
	{ 
	  printf("Executing %s...\n", argv[0]);
	  printf("**invalid command line flag >%c<>%c<\n", argv[argcount][0], argv[argcount][1]);
	  printf("Puta, puta, puta, try -h for help.\n");
	  exit(0);
	}
    }
}
