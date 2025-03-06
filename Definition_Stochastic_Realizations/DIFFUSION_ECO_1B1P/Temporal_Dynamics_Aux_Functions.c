#include <MODEL.h>

extern gsl_rng * r;   /* Global generator variable defined at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

void assert_Total_Population (Parameter_Table * Table, double * Y);
void Print_Discrete_Probability_Distribution(Parameter_Table * Table, int Event, int x);

#define ASSERTION_TRUE

void Carrying_Capacity_Control( int Event, Parameter_Table * Table, int x, int jS, double Y, int J )
{
  int i, Q, nS;
  int Upper_Boundary;
  treenode * nodeDum; 

  Community ** Patch = Table->Patch_System;
  double K_R         = (double)Table-> K_R; 

  #if defined ASSERTION_TRUE

    nS = jS%Table->LOCAL_STATE_VARIABLES;   /*Strain ID of the affected species */

    assert( nS == 0 || nS ==1 );
  

    Upper_Boundary = 0;                     /* If K_R is overpassed, the upper 
                                               limit has been overcome 
                                            */
    if ( Y >= K_R )                    Upper_Boundary = 1;
    if ( J >=  (int)K_R)               Upper_Boundary = 1;
    if ( Patch[x]->n[nS] >= (int)K_R ) Upper_Boundary = 1;                                           

    if( Upper_Boundary == 1 ) {
      printf (" Carrying Capacity control is not passed...\n");
      printf (" Warning: some population is over maximum capacity of the local community, but it should not!!!\n");
      printf (" Event No %d:", Event); printf("... in patch No %d\n", x);
      printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
      printf ("J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
      printf ("n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);
      for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
        printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
      for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++)
        printf ("Variable: %d\t Population: %d\n", i, Patch[x]->n[i]);

      if( Table->No_of_CELLS <= 4)   
        Print_Meta_Community_Patch_System (Table);

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      if(Table->No_of_CELLS <= 4) 
        nodeDum = leafPrint(Table->Treeroot); /* Printing leaf values */
        
      printf("Time: %g\t Next_Time(Event=%d) = %g\n", 
              Table->T->Rate->Stochastic_Time, Event, Table->Treeroot->value);
        // Print the structure of the priority queue tree for debugging purposes
        // printtree(Table->Treeroot);       
      #endif 

      /* If positivity control fails, the programs exits!!! */
      exit(0);
    }

  #endif
}

void Positivity_Control( int Event, Parameter_Table * Table,
			                   int x, int jS, double Y, int J )
{
  int i, Q, nS, Non_Positive;
  Community ** Patch = Table->Patch_System;
  treenode * nodeDum;

#if defined ASSERTION_TRUE

  nS = jS%Table->LOCAL_STATE_VARIABLES;    /*Strain ID of the affected species */

  assert( nS == 0 || nS == 1 );

  Non_Positive = 0;

  if ( Y <= 0.0 )            Non_Positive = 1;
  if ( J <= 0 )              Non_Positive = 1;
  if ( Patch[x]->n[nS] <= 0) Non_Positive = 1;

  if( Non_Positive == 1 ) {
    printf (" Positivity control is not passed...\n");
    // if( Table->No_of_CELLS <= 4)   
    Print_Meta_Community_Patch_System (Table);

    printf (" Event No %d:", Event); printf("... in patch No %d\n", x);
    printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
    printf (" J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
    printf (" n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);
    Print_Discrete_Probability_Distribution(Table, Event, x);
    printf (" Warning: some population is negative or zero, but it should not!!!\n");
      
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      if(Table->No_of_CELLS <= 4) 
        nodeDum = leafPrint(Table->Treeroot); /* Printing leaf values */
      
      printf("Time: %g\t Next_Time(Event=%d) = %g\n", 
              Table->T->Rate->Stochastic_Time, Event, Table->Treeroot->value);
      // printtree(Table->Treeroot);       
    #endif 
    
    /* If positivity control fails, the programs exits!!! */
    exit(0);
  }

#endif
}

int Carrying_Capacity_Control_Check( int Event, Parameter_Table * Table, int x, int jS, double Y, int J )
{
  treenode * nodeDum;
  int success; 
  int i, Q, nS;
  int Upper_Boundary;

  Community ** Patch = Table->Patch_System;
  double K_R         = (double)Table-> K_R; 

  success = 1;  

    nS = jS%Table->LOCAL_STATE_VARIABLES;   /*Strain ID of the affected species */

    assert( nS == 0 || nS ==1 );
  

    Upper_Boundary = 0;                     /* If K_R is overpassed, the upper 
                                               limit has been overcome 
                                            */
    if ( Y >= K_R )                    Upper_Boundary = 1;
    if ( J >=  (int)K_R)               Upper_Boundary = 1;
    if ( Patch[x]->n[nS] >= (int)K_R ) Upper_Boundary = 1;                                           

    if( Upper_Boundary == 1 ) {
      printf (" Carrying Capacity control is not passed...\n");
      printf (" Warning: some population is over maximum capacity of the local community, but it should not!!!\n");
      printf (" Event No %d:", Event); printf("... in patch No %d\n", x);
      printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
      printf ("J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
      printf ("n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);
      for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
        printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
      for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++)
        printf ("Variable: %d\t Population: %d\n", i, Patch[x]->n[i]);

      if( Table->No_of_CELLS <= 4)   
        Print_Meta_Community_Patch_System (Table);

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      if(Table->No_of_CELLS <= 4) 
        nodeDum = leafPrint(Table->Treeroot); /* Printing leaf values */
      
      printf("Time: %g\t Next_Time(Event=%d) = %g\n", 
              Table->T->Rate->Stochastic_Time, Event, Table->Treeroot->value);
        // Print the structure of the priority queue tree for debugging purposes
        // printtree(Table->Treeroot);       
      #endif 
 
      success = 0;
    }

  return(success);
}


void Print_Discrete_Probability_Distribution(Parameter_Table * Table, int Event, int x)
{
  /* Discrete probability distribution at the local patch */
  int i;
  Community ** Patch = Table->Patch_System;

  for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
      printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
}

