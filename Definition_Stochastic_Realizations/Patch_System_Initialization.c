#include <MODEL.h>

void Patch_System_Initialization (Community ** PATCH, Parameter_Table * Table, double * y_INI)
{
  int i,j, S;
  double x;
  
  /* Populations are initialized in agreement with y_INI  */

  /* S is the number of variables required to define the state of a single patch */
  S = Table->No_of_SPECIES; 

  x = 0.0; 
  for (j = 0; j < Table->No_of_CELLS; j++) {
    for(i=0; i < S; i++) {
      PATCH[j]->n[i] = (int)y_INI[i + j*S];
      x += y_INI[i + j*S];
    }
  }
                                            
  printf("Initial Total Population (per Species): %g\t", Table->INITIAL_TOTAL_POPULATION);

  printf("Total Community Size  (across Species): %g\n", x);

  getchar(); 
}
