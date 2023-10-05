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
  //Print_Press_Key(1,0,".");
  
  free(a);
  return(p-1);
}

int Choose_Village_Binary_Tree(double max_Probability, Community ** Pop, Parameter_Model * Par)
{
  /* This function depends on a function from treenode library called choose_Indvidual_Event() 
     (see Treenode.c)
  */
  int i, p;
  int No_of_Villages;

  No_of_Villages = Par->No_of_CELLS;

  double x = Par->Treeroot->value * RANDOM;
  p = choose_Individual_Event(Par->Treeroot, x);  /* From treenode.c library */
  
  /* p should go from (0) to (No of Villages-1) */
  assert(p >= 0 && p < No_of_Villages); 
  /* The total rate of change for each patch is correctly stored at leave level */
  assert(Pop[p]->ratePatch == Par->Leaves[p]->value);
  
  return(p);
}

void Choose_Village_and_Event_Binary_Tree( double max_Probability, Community ** Pop,
                                           Parameter_Model * Par, 
                                           int * patch, int * event )
{
  /* This function depends on a function from treenode library called choose_Indvidual_Event() 
     (see Treenode.c)
  */
  int i, p, n;
  int No_of_EVENTS, No_of_TREE_LEVELS;
  double S; 
 
  No_of_EVENTS  = Par->TOTAL_No_of_EVENTS;
  No_of_TREE_LEVELS = Par->No_of_TREE_LEVELS; 
 
  double x = Par->Treeroot->value * RANDOM;
  p = choose_Individual_Event(Par->Treeroot, x);  /* From treenode.c library */
  
  /* p should go from (0) to (TOTAL_GRAND_No_of_EVENTS-1) */
  assert(p >= 0 && p < Par->TOTAL_GRAND_No_of_EVENTS);

  * event = p%No_of_EVENTS; 
  * patch = p/No_of_EVENTS; 

  /* The rate of each event is correctly stored at leave level */
  assert(Par->Leaves[p]->value == Pop[*patch]->rToI[*event]);

  #if defined REUSE_RANDOM_NUMBER 
    /* Possibility to store a generated random number for further re-use!!! */
    n = No_of_TREE_LEVELS -1 ; 
    S = 0.0;
    partial_sums_upto_p(Par->Treeroot, p, n, &S);
    Par->Time->Rate->Reusable_Random_Number = (x - S)/Par->Leaves[p]->value; 
  #endif
}

void Choose_Village_and_Event_Priority_Queu(double max_Probability, Community ** Pop,
                                       Parameter_Model * Par, 
                                       int * patch, int * event )
{
  int i, p;
  Community * P;
  int No_of_EVENTS;

  No_of_EVENTS  = Par->TOTAL_No_of_EVENTS;

  p = Par->Treeroot->index; /* Just looking it up (at the root node, the heap of the 
                               priority queue */

  /* p should go from (0) to (No of Villages-1) */
  assert(p >= 0 && p < Par->TOTAL_GRAND_No_of_EVENTS);

  * event = p%No_of_EVENTS; 
  * patch = p/No_of_EVENTS;
}
                                        
