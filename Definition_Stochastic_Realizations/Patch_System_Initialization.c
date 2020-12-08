#include <MODEL.h>

void Patch_System_Initialization (Community ** PATCH, Parameter_Table * Table, double * y_INI)
{
  int i,j, Q; 
  
  /* Populations are initialized in agreement with y_INI  */

  /* Q is the number of variables required to define the state of a single patch */
  Q = Table->No_of_SPECIES; 

  for (j = 0; j < Table->No_of_CELLS; j++) {

    for(i=0; i<Q; i++) 
      PATCH[j]->n[i] = (int)y_INI[i + j*Q]; 
    
  }  
}
