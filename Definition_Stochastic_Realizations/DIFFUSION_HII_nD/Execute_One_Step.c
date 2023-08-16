/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                             David Alonso, 2010 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

extern gsl_rng * r;   /* Global generator (define at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

void assert_Total_Population (Parameter_Table * Table, double * Y);

#define ASSERTION_TRUE

void Execute_One_Step(Community ** SP,
		      Parameter_Table * Table,
		      double max_Probability,
		      int * Event, int * x_Patch)
{
  int x, y, n_Event, k, n, n_Event_Sp, Sp, j;
  Community * Patch;
  Parameter_Model * P = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  /* Hierarchic procedure to find the even to occur... */
  /* The event occurs in one of the local populations  */
  if(P->No_of_CELLS == 1) x = y = 0;
  else                    x = y = Choose_Village(max_Probability, SP, P);

  /* x and y will end up differing only in case there is movemnt event!!! */

  Patch = SP[x];  /* x represents the chosen patch undegoing a change. */

  if(Table->TOTAL_No_of_EVENTS > 1) {
    n_Event = Discrete_Sampling(Patch->rToI, Table->TOTAL_No_of_EVENTS) - 1; 
    /* n_Event should be bounded 
       between 0 and Table->No_of_RESOURCES*Table->No_of_EVENTS -1 
    */
    n_Event_Sp = n_Event%Table->No_of_EVENTS;
    Sp         = n_Event/Table->No_of_EVENTS;
  }
  else {
    printf(" The total number of events that potentially could happen patch %d\n", x);
    printf(" is zero??? (TOTAL_No_of_EVENTS = %d)\n", Table->TOTAL_No_of_EVENTS);
    printf(" Something very wrong with your code\n");
    printf(" The program will exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  }

  j = x*Table->No_of_RESOURCES + Sp;

  assert( P->No_of_CELLS == 1);
  assert(n_Event < Table->TOTAL_No_of_EVENTS-2); 
  /* No movement between patches or immgration from outsite the system is allowed */
  // Print_Meta_Community_Patch_System (Table);
   
  if (n_Event < Table->TOTAL_No_of_EVENTS-2) { 
    
    switch( n_Event_Sp )
    {
    case 0: /* Consumer Consumption of Resource and dimmer formation */ /* R + A ---> RA */
      
      Y[j]++; J[j]++;  Patch->n[Sp]++;            /* Dimer formation */

      Table->TOTAL_No_of_FREE_CONSUMERS--;
      
      break;

    case 1: /* Dimmer degratation back into a free individual consumers */ /* RA ---> A  */
      
      Positivity_Control( n_Event, Table, x, j, Y[j], J[j] );
      Y[j]--; J[j]--;  Patch->n[Sp]--;

      Table->TOTAL_No_of_FREE_CONSUMERS++;

      break;

    default:
    /* Something is very very wrong!!! */
      printf("The number of event occurring should be 0 and %d*%d-1\n", 
        Table->No_of_RESOURCES, Table->No_of_EVENTS);
      printf("Event to Occur = %d\n", n_Event);
      Print_Press_Key(1,0,".");
      exit(0);
    }
  }
  else {
    
    /* DIFFUSION_HII_nD will never go through this else, because 
       no movement between patches or immgration from outsite the system is allowed 
    */  
    k = n_Event - Table->No_of_RESOURCES*Table->No_of_EVENTS;

    switch( k )
    {
    case  0:  /* Out-Migration (A --> A-1) from patch x and some other patch gains one */
              /* Out A: only diffusion of free consumers    */
              /* 'A' is not defined properly for this model */
              /* 'A' should be equal to 
                Table->No_of_RESOURCES*Table->No_of_EVENTS, 
                the index of the variable representing free consumers 
              */  
      Positivity_Control( n_Event, Table, x, A, Y[A], J[A] );
      Y[A]--; J[A]--; Patch->n[Table->A]--;

      y = Some_Other_Patch_Population_Increase(x, Table->A, Table);

      break;

    case  1:  /* External Immigration from outside the system (A --> A+1) *//* External Imm A */

      Y[A]++; J[A]++;  Patch->n[Table->A]++;

      break;

    default:
    /* Something is very very wrong!!! */
      printf("At this point, the number of event occurring should be eitenr 0 or 1\n");
      printf("Event to Occur = %d\n", n_Event);
      Print_Press_Key(1,0,".");
      exit(0);
    }
  }

  (*Event) = n_Event;  x_Patch[0] = x;  x_Patch[1] = y;
}

int Some_Other_Patch_Population_Increase(int x, int Sp, Parameter_Table * Table)
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

  j = Sp + Patch[x]->Patch_Connections[n_Patch]*Table->LOCAL_STATE_VARIABLES;

  Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;

  return(Patch[x]->Patch_Connections[n_Patch]);
}

void Positivity_Control( int Event, Parameter_Table * Table, int x, int jS, double Y, int J)
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
    printf (" Event No %d:", Event);
    printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
    printf ("J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
    printf ("n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);

    for(i=0; i<Table->TOTAL_No_of_EVENTS; i++) 
      printf ("Event: %d\t Rate of Event No %d: %g\n", Event, i, Patch[x]->rToI[i]);
    
    for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++) 
      printf ("Varible: %d\t Population: %d\n", i, Patch[x]->n[i]);

    Print_Meta_Community_Patch_System (Table);
    // exit(0);
  }

#endif
}
