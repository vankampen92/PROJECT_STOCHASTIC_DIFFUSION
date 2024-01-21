#include <MODEL.h>

void Initial_Condition_Boundary_Values(int j, double *value_0, double *value_1,
                                Parameter_Table * Table)
{
/* Definition of the state vector numerical order, from 0 to K, of model variables */
   #include <Model_Variables_Code.Include.c>
  
   if ( j < Table->MODEL_STATE_VARIABLES ) {
     (*value_1) = Table->INITIAL_TOTAL_POPULATION; /* Total Initial Population per Sp */
     (*value_0) = 0.0;       /* 0 */
   }
   else{
     printf(".... INVALID VARIABLE KEY [key = %d]\n", j);
     printf(".... The permited correspondences are:\n");
     printf(".... The program will exit\n");
     exit(0);
   }
}
