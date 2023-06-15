#include <MODEL.h>

extern gsl_rng * r;   /* Global generator (define at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

int Choose_Village(double max_Probability, Community ** Pop, Parameter_Model * Par)
{
  /* This function depends on Discrete_Sampling() */
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

int Choose_Village_Binary_Tree(double max_Probability, Community ** Pop, Parameter_Model * Par)
{
  /* This function depends on a function from treenode library called choose_Indvidual_Event() */
  int i, p;
  Community * P;
  int No_of_Villages;

  No_of_Villages = Par->No_of_CELLS;

  double x = Par->Treeroot->value * RANDOM; 
  p = choose_Individual_Event(Par->Treeroot, x);  /* From treenode.c library */
  
  /* p should go from (0) to (No of Villages-1) */
  assert(p >= 0 && p < No_of_Villages); 
  
  return(p);
}