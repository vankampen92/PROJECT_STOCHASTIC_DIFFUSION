/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                             David Alonso, 2010 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

extern gsl_rng * r;   /* Global generator variable defined at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

void assert_Total_Population (Parameter_Table * Table, double * Y);
void Print_Discrete_Probability_Distribution(Parameter_Table * Table, int Event, int x);

#define ASSERTION_TRUE

#ifndef BINARY_TREE_SUPER_OPTIMIZATION
#ifndef PRIORITY_QUEU_SUPER_OPTIMIZATION
  #define EVENT_DISCRETE_SAMPLING 
#endif
#endif

void Execute_One_Step(Community ** SP,
		                  Parameter_Table * Table,
		                  double max_Probability,
		                  int * Event, int * x_Patch)
{
  int x, y, n_Event, Index_Patch, Index_Event, n, n_Event_Sp, Sp, j;
  Community * Patch;
  Parameter_Model * P = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  /* Hierarchic procedure to find the even to occur... */
  /* The event occurs in one of the local populations  */
  if(P->No_of_CELLS == 1) {
    x = y = 0;
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Choose_Village_and_Event_Binary_Tree(max_Probability, SP, P, 
                                           &Index_Patch, &Index_Event);
      x = y   = Index_Patch; 
      n_Event = Index_Event;
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Choose_Village_and_Event_Priority_Queu(max_Probability, SP, P, 
                                             &Index_Patch, &Index_Event);
      x = y   = Index_Patch; 
      n_Event = Index_Event;

      assert(Index_Event == Table->Treeroot->index );
    #endif
    assert(x == 0); 
    assert(y == 0);
  }
  else {
    #if defined BINARY_TREE_OPTIMIZATION
        x = y = Choose_Village_Binary_Tree(max_Probability, SP, P);
    #elif defined BINARY_TREE_SUPER_OPTIMIZATION
        Choose_Village_and_Event_Binary_Tree(max_Probability, SP, P, 
                                             &Index_Patch, &Index_Event);
        x = y   = Index_Patch; 
        n_Event = Index_Event;
    #elif defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Choose_Village_and_Event_Priority_Queu(max_Probability, SP, P, 
                                               &Index_Patch, &Index_Event);
        x = y   = Index_Patch; 
        n_Event = Index_Event;

        assert( (Index_Patch * Table->TOTAL_No_of_EVENTS + Index_Event) == Table->Treeroot->index );

        assert( Table->Tree_Node_Index[Table->Treeroot->index] == Table->Treeroot );
    #else
        x = y = Choose_Village(max_Probability, SP, P);
    #endif
  }
  /* When this function 'returns', 
     x and y will end up being different only in case there is movemnt event!!! 
  */
  Patch = SP[x];  /* x represents the chosen patch undegoing a change. */
  #if defined EVENT_DISCRETE_SAMPLING
    if(Table->TOTAL_No_of_EVENTS > 1) 
      n_Event = Discrete_Sampling(Patch->rToI, Table->TOTAL_No_of_EVENTS) - 1; 
      /* 0, ..., Tablel->TOTAL_No_of_EVENTS-1 */
  #endif
  
  if (Table->TOTAL_No_of_EVENTS == 0) {
    printf(" The potential number of events that potentially could happen in patch %d\n", x);
    printf(" is zero??? (TOTAL_No_of_EVENTS = %d)\n", Table->TOTAL_No_of_EVENTS);
    printf(" If it is, there is something very wrong with your code\n");
    printf(" The program will exit!!!\n");
    Print_Press_Key(1,1,"Error in Execute_One_Step() funciton");
    exit(0);
  }

  // Print_Discrete_Probability_Distribution(Table, n_Event, x); /* Comment out if it works */

  W   = x*Table->LOCAL_STATE_VARIABLES + Table->W;
  Q    = x*Table->LOCAL_STATE_VARIABLES + Table->Q;
  F    = x*Table->LOCAL_STATE_VARIABLES + Table->F;
  WF   = x*Table->LOCAL_STATE_VARIABLES + Table->WF;

  assert( n_Event < Table->TOTAL_No_of_EVENTS );
  assert( x       < Table->No_of_CELLS );

  // Print_Meta_Community_Patch_System (Table);

  switch( n_Event )
    {
    case  0:  /* Out-Migration (W --> W - 1) from patch x to some other patch */         /* Outmigraiton (W) */
      Positivity_Control( 0, Table, x, W, Y[W], J[W] );

      Y[W]--; J[W]--; Patch->n[Table->W]--;

      y = Some_Other_Patch_Population_Increase(x, Table->W, Table);

      break;

    case  1:  /* External Immigration from outside the system (W --> W+1) */             /* External Imm W */

      Y[W]++; J[W]++;  Patch->n[Table->W]++;

      break;

    case  2:  /* Local Death of Workers (W --> W-1)  */     /* Death (W) */
      Positivity_Control( 2, Table, x, W, Y[W], J[W] );

      Y[W]--; J[W]--;  Patch->n[Table->W]--;

      break;

    case  3:  /* Worker Production by Queens (W --> W+1) */                                     /* Production (W) */

      Y[W]++; J[W]++;  Patch->n[Table->W]++;

      break;

    case  4:  /* Nest Establishment (Q ---> Q+1 and W ---> W-1)  */ /* Establishment (W) */
      Positivity_Control( 4, Table, x, W, Y[W], J[W] );

      Y[W]--; J[W]--;  Patch->n[Table->W]--;
      Y[Q]++;  J[Q]++;   Patch->n[Table->Q]++;

      break;

    case  5:  /* Death of Queens (Q ---> Q-1)  */                                       /* Death (Q) */
      Positivity_Control( 5, Table, x, Q, Y[Q], J[Q] );

      Y[Q]--; J[Q]--;  Patch->n[Table->Q]--;

      break;

    case  6:  /* Out-Migration (F --> F-1) from patch x and some other patch 
                gains one */  /* Out of Flies */
      Positivity_Control( 6, Table, x, F, Y[F], J[F] );

      Y[F]--; J[F]--; Patch->n[Table->F]--;

      y = Some_Other_Patch_Population_Increase(x, Table->F, Table);

      break;

    case  7:  /* External Immigration from outside the system (F --> F+1) */               /* External Imm (F) */

      Y[F]++; J[F]++;  Patch->n[Table->F]++;

      break;

    case  8:  /* Local Death  (F --> F-1)  */                                               /* Death (F) */
      Positivity_Control( 8, Table, x, F, Y[F], J[F] );

      Y[F]--; J[F]--;  Patch->n[Table->F]--;

      break;

    case 9: /* Attack: Fly attack on workers */ /* W + F --->  WF + F */
      Positivity_Control( 9, Table, x, W, Y[W], J[W] );

      Y[W]--; J[W]--; Patch->n[Table->W]--;
      
      Y[WF]++; J[WF]++;  Patch->n[Table->WF]++;
      break;

    case 10: /* Fly Laval Development into Adult Flies */  /* WF ---> F */
      Positivity_Control( 10, Table, x, WF, Y[WF], J[WF] );

      Y[WF]--; J[WF]--;  Patch->n[Table->WF]--;
      Y[F]++;  J[F]++;   Patch->n[Table->F]++;

      break;
      
    case 11: /* Local Death of Parasitized Workers (WF)   /* WF  ---> WF - 1  */
      Positivity_Control( 11, Table, x, WF, Y[WF], J[WF] );
      
      Y[WF]--; J[WF]--;  Patch->n[Table->WF]--;

      break;
      
    default:
    /* Something is very very wrong!!! */
      printf("The number of event occurring should be between 0 and 12\n");
      printf("Event to Occur = %d\n", n_Event);
      Print_Press_Key(1,0,".");
      exit(0);
    }

  (*Event) = n_Event;  x_Patch[0] = x;  x_Patch[1] = y;
}

int Some_Other_Patch_Population_Increase(int x, int Sp,
					                               Parameter_Table * Table)
{
  /* Input:

     .  x: Patch label of the patch from which the individuals goes out
     . Sp: Species;
     . Table: Parameter Table

     Output:

     . y: Patch label of the patch receiving the immigrant

  */
  int k, j, n_Patch;
  Community ** Patch = Table->Patch_System;

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  n_Patch = Discrete_Sampling(Patch[x]->Out_Migration_Vector[Sp], Patch[x]->No_NEI) - 1;

  j = Sp + Patch[x]->Patch_Connections[n_Patch]*Table->LOCAL_STATE_VAQIABLES;

  Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;

  return(Patch[x]->Patch_Connections[n_Patch]);
}

void Positivity_Control( int Event, Parameter_Table * Table,
			                   int x, int jS, double Y, int J)
{
  int i, Q, nS, Non_Positive;
  Community ** Patch = Table->Patch_System;

#if defined ASSEQTION_TRUE

  nS = jS%Table->LOCAL_STATE_VARIABLES;

  Non_Positive = 0;

  if ( Y <= 0.0 )            Non_Positive = 1;

  if ( J <= 0 )              Non_Positive = 1;

  if ( Patch[x]->n[nS] <= 0) Non_Positive = 1;

  if( Non_Positive == 1 ) {
    printf (" Positivity control is not passed...\n");
    printf (" Warning: some population is negative or zero, but it should not!!!\n");
    printf (" Event No %d:", Event); printf("... in patch No %d\n", x);
    printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
    printf ("J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
    printf ("n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);
    for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
      printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
    for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++)
      printf ("Varible: %d\t Population: %d\n", i, Patch[x]->n[i]);

    Print_Meta_Community_Patch_System (Table);

    treenode * nodeDum = leafPrint(Table->Treeroot);

    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      printf("Time: %g\t Next_Time(Event=%d) = %g\n", 
              Table->T->Rate->Stochastic_Time, Event, Table->Treeroot->value);
      printtree(Table->Treeroot);       
    #endif 

    exit(0);
  }

#endif
}

void Print_Discrete_Probability_Distribution(Parameter_Table * Table, int Event, int x)
{
  int i;
  Community ** Patch = Table->Patch_System;

  for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
      printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
}