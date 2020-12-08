#include <MODEL.h>

int Choose_Village(double max_Probability, Community ** Pop, Parameter_Model * Par)
{
  int i,p;
  Community * P;
  int No_of_Villages;

  No_of_Villages = Par->No_of_CELLS;

  double * a = (double *)malloc(No_of_Villages * sizeof(double) );
  
  for(i=0; i < No_of_Villages; i++) {
    P  = Pop[i];
    a[i] = P->ratePatch;
  }

  p = Discrete_Sampling(a, No_of_Villages);
  //p = Discrete_Sampling_Rejection_Method(max_Probability, a, No_of_Villages);

  //for(i=0; i < No_of_Villages; i++) printf("a[%d] = %g----", i, a[i]);
  //Press_Key();
  
  free(a);
  return(p-1);
}
