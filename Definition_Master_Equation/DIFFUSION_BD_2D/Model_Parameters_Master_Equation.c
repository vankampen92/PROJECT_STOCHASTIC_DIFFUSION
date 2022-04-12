#include <MODEL.h>

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				      int * No_of_CONFIGURATIONAL_STATES,
				      int * n_DIMENSION,
				      int * n_x, int * n_y, int * n_z)
{
  
  
  * n_DIMENSION = 2;
  * n_z         = 1;
  * n_y         = Table->TOTAL_No_of_CONSUMERS;
  * n_x         = FEEDING_EVENT_FACTOR * (int)N; 
  
  * No_of_CONFIGURATIONAL_STATES = Table->TOTAL_No_of_CONSUMERS * ( * n_x ) ;
}










