#include <MODEL.h>

void Initial_Condition_Master_Equation( Parameter_Table * Table, double * y_INI )
{
  /* All Consumers are free at time 0. No consumers handling resources at the 
     initial time.
  */
  int i, i_0, No_of_CONFIGURATIONAL_STATES; 
  
  /* Vector of Handling Consumers */
  int * n = (int *)calloc(Table->No_of_RESOURCES, sizeof(int));

  No_of_CONFIGURATIONAL_STATES = Table->MEq->No_of_CONFIGURATIONAL_STATES;
  Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = Table->TOTAL_No_of_CONSUMERS;
  Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = 0.0;

  /* Notice that we assign probability 1 to the configurational state characterized 
     by the absence of handlig consumers (on any resource):
     
                        P(n, t) = 1, where n = (0, ..., 0)
  */
     
  i_0 = Configuration_to_i_Map(n, Table->No_of_RESOURCES, Table->TOTAL_No_of_CONSUMERS);

  assert(i_0 == 0); /* by definition */
  
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {
    
    if( i == i_0 ) y_INI[i] = 1.0;
    else           y_INI[i] = 0.0;       
    
  }

  free(n);
}

							 
  
  
  
