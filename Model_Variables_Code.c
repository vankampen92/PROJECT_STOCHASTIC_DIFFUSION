#include <MODEL.h>

void Model_Variables_Code_into_Parameter_Table (Parameter_Table * Table)
{
  int i, j, n;
  int TYPE_of_MODEL;

  TYPE_of_MODEL = Table->TYPE_of_MODEL;

  switch( TYPE_of_MODEL )
    {

    case 0: /* DIFFUSION * * * * * * * * * * * * * * * * * * * * * * */

      /* Number of events that can occur to a single Species: */
      Table->No_of_EVENTS       = 1;  /* (Only Diffusion)         */
      Table->TOTAL_No_of_EVENTS = Table->No_of_EVENTS * Table->No_of_RESOURCES;
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->No_of_RESOURCES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */

      n = 0; 
      Table->Index[n++] = 0; /* Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */  

      Table->TOTAL_No_of_MODEL_PARAMETERS = n; 
      
      break;

    case 1: /* DIFFUSION_S_RESOURCES * * * * * * * * * * * * * * * * * * * * * * */

      /* Number of events that can occur to a single Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = Table->No_of_EVENTS * Table->No_of_RESOURCES;
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->No_of_RESOURCES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */

      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      if ( Table->No_of_RESOURCES > 0 ) {
	Table->Index[n++] = 6;
	Table->Index[n++] = 7; 
      }
      if ( Table->No_of_RESOURCES > 1 ) {
	Table->Index[n++] = 8;
	Table->Index[n++] = 9; 
      }
      Table->Index[n++]   = 10; /* No_of_RESOURCES */ 
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      
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
