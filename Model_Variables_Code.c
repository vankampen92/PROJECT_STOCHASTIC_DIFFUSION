#include <MODEL.h>

void Model_Variables_Code_into_Parameter_Table (Parameter_Table * P)
{
  int i, j, n;
  int TYPE_of_MODEL;

  TYPE_of_MODEL = P->TYPE_of_MODEL;

  switch( TYPE_of_MODEL )
    {

    case 0: /* DIFFUSION * * * * * * * * * * * * * * * * * * * * * * */

      /* Number of events that can occur to a single Species: */
      P->No_of_EVENTS       = 1;  /* (Only Diffusion)         */
      P->TOTAL_No_of_EVENTS = P->No_of_EVENTS * P->No_of_SPECIES;
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	       for(j=0; j<P->No_of_SPECIES; j++)
	        //P->n[i] = n++;
	         n++;
	    
      /* Conventions */
      P->K   = n-1;     /* Label last class            */

      break;

    default:
      printf(" This TYPE_of_MODEL (%d) code is not defined.\n", TYPE_of_MODEL);
      printf(" Check input argument list\n");
      exit(0);
   }
  /* Conventionally, the last label in the argument list of

     Model_Variables_Code_into_Parameter_Table (...),

     (*K), should be the label of the last model state variable.
     Then ( * K) + 1 becomes de total number of dynamic variables.
  */
}
