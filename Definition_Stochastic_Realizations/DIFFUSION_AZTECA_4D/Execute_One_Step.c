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

// #define ASSERTION_TRUE

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

      y = Other_Patch_Population_Increase(0, x, W, Table->W, Table);
      /* If immigration fails, the loss of an invidual in the source patch 
         is equivalent to a local death event. In fact, the individual dies 
         while trying to settle in the new patch */
      if (y == x) n_Event = 2; 

      break;

    case  1:  /* External Immigration from outside the system (W --> W+1) */             /* External Imm W */

      assert(Table->Lambda_R_0 > 0.0); 
      
      Boundness_Control_Local_Patch( 1, Table, x, W, Y[W], J[W] );
      Y[W]++; J[W]++;  Patch->n[Table->W]++;

      break;

    case  2:  /* Local Death of Workers (W --> W-1)  */     /* Death (W) */
      Positivity_Control( 2, Table, x, W, Y[W], J[W] );
      Y[W]--; J[W]--;  Patch->n[Table->W]--;

      break;

    case  3:  /* Worker Production by Queens (W --> W+1) */                                     /* Production (W) */
      Boundness_Control_Local_Patch( 3, Table, x, W, Y[W], J[W] );
      Y[W]++; J[W]++;  Patch->n[Table->W]++;

      break;

    case  4:  /* Nest Establishment (Q ---> Q+1 and W ---> W-1)  */ /* Establishment (W) */
      Positivity_Control( 4, Table, x, W, Y[W], J[W] );
      Y[W]--; J[W]--;  Patch->n[Table->W]--;

      Boundness_Control_Local_Patch( 4, Table, x, Q, Y[Q], J[Q] );
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

      y = Other_Patch_Population_Increase(6, x, F, Table->F, Table);

      /* If immigration fails, the loss of an individual in the source patch 
         is equivalent to a local death event. In fact, the individual 'dies' 
         while trying to settle in the new patch, because the new patch is 
         full (at carrying capacity) */
      if (y == x) n_Event = 8;

      break;

    case  7:  /* External Immigration from outside the system (F --> F+1) */               /* External Imm (F) */

      assert(Table->Lambda_C_0 > 0.0); 
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

int Other_Patch_Population_Increase(int Event, int x, int n, int Sp,
					                          Parameter_Table * Table)
{
  /* Input:
     . Event, process responsable for that possible increase 
     . x: Patch label of the patch from which the individuals goes out
     . n:  Index identifying the patch from which the individual moves 
           and the species in the arrays Vector_Model_Variables[] and 
           Vector_Model_Variables_Int[] of the main Table structure. 
     . Sp: Species label of the individual that jumps.
     . Table: Parameter Table

     Output:
     . y: Patch label of the patch receiving the immigrant

  */
  int k, j, y, n_Patch;
  Community ** Patch = Table->Patch_System;

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  n_Patch = Discrete_Sampling(Patch[x]->Out_Migration_Vector[Sp], Patch[x]->No_NEI) - 1;

  y = Patch[x]->Patch_Connections[n_Patch];

  j = Sp + Patch[x]->Patch_Connections[n_Patch]*Table->LOCAL_STATE_VARIABLES;

  /* n_W (number of  workers) has to be always lower or equal to K_W */
  if( Sp == 0 || Sp == 1) {
    y = Boundness_Control( Event, Table, x, y, Sp, j, Y[j], J[j] );
    if (y == x) { 
      /* No increase in the receiving patch is possible */
       y = x;  
    }
    else { 
      Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;
    }
  }
  else { 
    Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;
  }  

  return(y);
}

int Boundness_Control( int Event, Parameter_Table * Table,
			                 int x_out, int x_in, int Sp, int jS, double Y, int J)
{
  /* Input: 
      . Event
      . Table
      . x_in, patch that receives the 'jumping' individual  (sink patch)
      . x_out, patch that individuals comes from            (source patch)
      . Sp, species label of the species involved
      . jS
      . Y and J are the quantities to be checked. They correspond to
          J = Table->Vector_Model_Variables_Int[jS], and 
          Y = Table->Vector_Model_Variables[jS]
     Output:
      . y 
  */
  int y; 
  int i, Q, nS, Non_Bounded;
  Community ** Patch = Table->Patch_System;
  double K; 

  /* K_W: Total Carrying Capacity (Workers), 
          and Table->K_R is the max No of Worker per Nest */
  double K_W = (double)Table->K_R * (double)Table->Lambda_C_1;  /* -HK  -H7 */
  /* K_Q: Total Max No of Nests (per local patch) */ 
  double K_Q = (double)Table->Lambda_C_1;                       /* -H7      */

  if (Sp == 0)
    K = K_W; 
  else if (Sp == 1)
    K = K_Q; 
  else {
    printf("Sp %s (index=%d) to increase is not bounded by any carrying capacity\n", 
            Table->Model_Variable_Name[jS], Sp);
    printf("but Boundness_Control() is called as if it were. Some must be wrong!!!\n");
    printf("The program will exit\n");
    exit(0);
  }

  nS = jS%Table->LOCAL_STATE_VARIABLES;
  y = x_in; 

  Non_Bounded = 0;
  if ( Y >= K )                       Non_Bounded = 1;
  if ( Patch[x_in]->n[nS] >= (int)K ) Non_Bounded = 1;

  if( Non_Bounded == 1 ) {
    printf (" Boundness control is not passed...\n");
    printf (" Warning: some population is higher that carrying capacity, but it should not!!!\n");
    printf (" Event No %d:", Event); printf("(reveiving patch No %d)\n", x_in);
    printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
    printf (" J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
    printf (" n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x_out]->n[nS]);
    printf (" Events, rates, and populations in source patch (%d):\n", x_out);
    for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
      printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x_out]->rToI[i]);
    for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++)
      printf ("Varible: %d\t Population: %d\n", i, Patch[x_out]->n[i]);

    Print_Meta_Community_Patch_System (Table);

    #if defined BINARY_TREE_OPTIMIZATION
      treenode * nodeDum = leafPrint(Table->Treeroot);
      printf("Time: %g\t Total Rate = %g\n", 
              Table->T->Rate->Stochastic_Time, Table->Treeroot->value);
      printtree(Table->Treeroot);  
    #endif

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      treenode * nodeDum = leafPrint(Table->Treeroot);
      printf("Time: %g\t Total Rate = %g\n", 
              Table->T->Rate->Stochastic_Time, Table->Treeroot->value);
      printtree(Table->Treeroot);  
    #endif

    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      treenode * nodeDum = leafPrint(Table->Treeroot);

      printf("Time: %g\t Next_Time(Event=%d) = %g\n", 
              Table->T->Rate->Stochastic_Time, Event, Table->Treeroot->value);
      printtree(Table->Treeroot);       
    #endif 

    /* The increase of the population on the receving patch x_in should fail... */
    /* The output of the function, y, will no longer be x_in, but the same patch
       the individual tried to leave */
    y = x_out;     
  }

  return(y);
}

void Boundness_Control_Local_Patch(int Event, Parameter_Table * Table, 
                                   int x, int jS, double Y, int J)
{
  int i, Q, nS, Non_Bounded;
  Community ** Patch = Table->Patch_System;
  double K; 
   
  #if defined ASSERTION_TRUE 
    nS = jS%Table->LOCAL_STATE_VARIABLES; /* Species being controled */
    Non_Bounded = 0;

    /* K_W: Total Carrying Capacity (Workers), 
          and Table->K_R is the max No of Worker per Nest */
    double K_W = (double)Table->K_R * (double)Table->Lambda_C_1;  /* -HK  -H7 */
    /* K_Q: Total Max No of Nests (per local patch) */ 
    double K_Q = (double)Table->Lambda_C_1;                       /* -H7      */

    if (nS == 0)      K = K_W; 
    else if (nS == 1) K = K_Q; 
    else {
      printf("Sp %s (index=%d) to increase is not bounded by any carrying capacity\n", 
            Table->Model_Variable_Name[jS], nS);
      printf("but Boundness_Control() is called as if it were. Some must be wrong!!!\n");
      printf("The program will exit\n");
      exit(0);
    }

    if ( Y >= K ) Non_Bounded = 1;
    if ( Patch[x]->n[nS] >= (int)K ) Non_Bounded = 1;

    if( Non_Bounded == 1 ) {
      printf (" Boundness control is not passed...\n");
      printf (" Warning: some population is larger than carrying capacity, but it should not!!!\n");
      printf (" Event No %d:", Event); printf("... in patch No %d\n", x);
      printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
      printf ("J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
      printf ("n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);
      for(i=0; i<Table->TOTAL_No_of_EVENTS; i++)
        printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
      for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++)
        printf ("Varible: %d\t Population: %d\n", i, Patch[x]->n[i]);

      Print_Meta_Community_Patch_System (Table);

      #if defined BINARY_TREE_OPTIMIZATION
        treenode * nodeDum = leafPrint(Table->Treeroot);
        printf("Time: %g\t Total Rate = %g\n", 
              Table->T->Rate->Stochastic_Time, Table->Treeroot->value);
        printtree(Table->Treeroot);  
      #endif

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        treenode * nodeDum = leafPrint(Table->Treeroot);
        printf("Time: %g\t Total Rate = %g\n", 
              Table->T->Rate->Stochastic_Time, Table->Treeroot->value);
        printtree(Table->Treeroot);  
      #endif

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        printf("Time: %g\t Next_Time(Event=%d) = %g\n", 
              Table->T->Rate->Stochastic_Time, Event, Table->Treeroot->value);
        printtree(Table->Treeroot);       
      #endif 
      
      exit(0);
    }
  #endif 
}
                       
void Positivity_Control( int Event, Parameter_Table * Table,
			                   int x, int jS, double Y, int J)
{
  int i, Q, nS, Non_Positive;
  Community ** Patch = Table->Patch_System;

#if defined ASSERTION_TRUE

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

    #if defined BINARY_TREE_OPTIMIZATION
      treenode * nodeDum = leafPrint(Table->Treeroot);
      printf("Time: %g\t Total Rate = %g\n", 
              Table->T->Rate->Stochastic_Time, Table->Treeroot->value);
      printtree(Table->Treeroot);  
    #endif

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      treenode * nodeDum = leafPrint(Table->Treeroot);
      printf("Time: %g\t Total Rate = %g\n", 
              Table->T->Rate->Stochastic_Time, Table->Treeroot->value);
      printtree(Table->Treeroot);  
    #endif

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
