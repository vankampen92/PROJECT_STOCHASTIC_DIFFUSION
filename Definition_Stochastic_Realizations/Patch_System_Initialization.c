#include <MODEL.h>

void Patch_System_Initialization (Community ** PATCH, Parameter_Table * Table, double * y_INI)
{
  int i,j, S;
  double x, x_S;

  /* Populations are initialized in agreement with y_INI  */

  /* S is the number of variables required to define the state of a single patch */
  S = Table->LOCAL_STATE_VARIABLES;

  x = 0.0;
  for(i=0; i < S; i++) {

    x_S = 0; ;
    for (j = 0; j < Table->No_of_CELLS; j++) {

      PATCH[j]->n[i] = (int)y_INI[i + j*S];
      x_S += y_INI[i + j*S];
      x   += y_INI[i + j*S];

#if defined DIFFUSION_1R1C
      if (i == 0) printf("Initial Population[Species R] in patch %d: %d\n", j, PATCH[j]->n[i]);
      if (i == 1) printf("Initial Population[Species A] in patch %d: %d\n", j, PATCH[j]->n[i]);
      if (i == 2) printf("Initial Population[Species AR] in patch %d: %d\n", j, PATCH[j]->n[i]);
      if (i == 3) printf("Initial Population[Species ARA] in patch %d: %d\n", j, PATCH[j]->n[i]);
#endif
    }

#if defined DIFFUSION_1R1C
    if (i == 0) printf("Initial Population over all patches [Species R]: %g\n", x_S);
    if (i == 1) printf("Initial Population over all patches [Species A]: %g\n", x_S);
    if (i == 2) printf("Initial Population over all patches [Species RA]: %g\n", x_S);
    if (i == 3) printf("Initial Population over all patches [Species ARA]: %g\n", x_S);
#endif
  }

#if defined DIFFUSION_HII_nD
  int TOTAL_No_of_HANDLING_CONSUMERS = 0;
  assert(Table->No_of_CELLS == 1);
  for(i=0; i < S; i++) {
    printf("Initial Number of Handling Consumers (in the %d-Patch) on the %d-Resource Type: %d\n",
           0, i, PATCH[0]->n[i]);
    TOTAL_No_of_HANDLING_CONSUMERS += PATCH[0]->n[i];
  }
  Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = TOTAL_No_of_HANDLING_CONSUMERS;
  Table->TOTAL_No_of_HANDLING_CONSUMERS        = TOTAL_No_of_HANDLING_CONSUMERS;
  Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0     = Table->TOTAL_No_of_CONSUMERS - TOTAL_No_of_HANDLING_CONSUMERS;
  Table->TOTAL_No_of_FREE_CONSUMERS            = Table->TOTAL_No_of_CONSUMERS - TOTAL_No_of_HANDLING_CONSUMERS;
#else
  printf("Initial Total Population (per Species): %g\t", Table->INITIAL_TOTAL_POPULATION);
  printf("Total Community Size  (across Species): %g\n", x);
#endif

  printf("Patch system successfully initialized\n");
  Press_Key();

#if defined VERBOSE
    Press_Key();
#endif
}
