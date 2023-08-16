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
  int x, y, n_Event, n, n_Event_Sp, Sp, j;
  Community * Patch;
  Parameter_Model * P = Table->P;

   /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  /* Hierarchic procedure to find the even to occur... */
  /* The event occurs in one of the local populations  */
  if(P->No_of_CELLS == 1) x = 0;
  else                    x = Choose_Village(max_Probability, SP, P);

  Patch = SP[x];  /* x represents the chosen patch undegoing a change. */

  if(Table->TOTAL_No_of_EVENTS > 1){
    n_Event = Discrete_Sampling(Patch->rToI, Table->TOTAL_No_of_EVENTS) - 1; /* 0, ..., 99 */
    n_Event_Sp = n_Event/Table->No_of_RESOURCES;
    Sp         = n_Event%Table->No_of_RESOURCES;
  }
  else {
    n_Event    = 0;
    n_Event_Sp = 0;
    Sp         = 0;
  }

  j = x*Table->No_of_RESOURCES + Sp;

  assert( n_Event < Table->No_of_RESOURCES );
  assert( n_Event_Sp == 0);

  // Print_Meta_Community_Patch_System (Table);

  switch( n_Event_Sp )
    {
    case  0:  /* Out-Migration (A --> A-1) and some other patch gains one */       /* Out S */
      Positivity_Control( 0, Table, x, j, Y[j], J[j] );

      Y[j]--; J[j]--;  Patch->n[Sp]--;

      y = Some_Other_Patch_Population_Increase(x, Sp, Table);

      break;
    /* case  1:  /\* In-Migration (A --> A+1) and some other patch loses one *\/ /\* In S *\/ */
    /*   Y[j]++; J[j]++;  Patch->n[Sp]++;                                                     */

    /*   Some_Other_Patch_Population_Decrease(x, Sp, Table);                                  */

    /*   break;                                                                               */

    default:
    /* Something is very very wrong!!! */
      printf("The number of event occurring should be between 0 and 0\n");
      printf("Event to Occur = %d\n", n_Event);
      Print_Press_Key(1,0,".");
      exit(0);
    }

#if defined ASSERTION_TRUE
  assert_Total_Population (Table, Y);
#endif

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

  j = Sp + Patch[x]->Patch_Connections[n_Patch]*Table->No_of_RESOURCES;

  Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;

  return(Patch[x]->Patch_Connections[n_Patch]);
}

void Positivity_Control( int Event, Parameter_Table * Table,
			 int x, int jS, double Y, int J)
{
  int Q, nS, Non_Positive;
  Community ** Patch = Table->Patch_System;

#if defined ASSERTION_TRUE

  nS = jS%Table->No_of_RESOURCES;

  Non_Positive = 0;

  if ( Y <= 0.0 )            Non_Positive = 1;

  if ( J <= 0 )              Non_Positive = 1;

  if ( Patch[x]->n[nS] <= 0) Non_Positive = 1;

  if( Non_Positive == 1 ) {
    printf (" Event No %d:", Event);
    printf (" Y[%s] = %g\t", Table->Model_Variable_Name[jS], Y);
    printf ("J[%s] = %d\t",  Table->Model_Variable_Name[jS], J);
    printf ("n[%s] = %d\n",  Table->Model_Variable_Name[jS], Patch[x]->n[nS]);
    exit(0);
  }

#endif
}

void assert_Total_Population (Parameter_Table * Table, double * Y)
{
  double N, N_Time_0;

  N = Total_Population (Y, Table ); /* Total Population is, in fact, Total Community Size */

  N_Time_0 = Table->No_of_RESOURCES * Table->No_of_INDIVIDUALS;

  if(N != N_Time_0)
    printf( "Total Community Size: %g\t S(%d) * No_of_INDIVIDUALS(%g) = %g\n",
	    N, Table->No_of_RESOURCES, Table->INITIAL_TOTAL_POPULATION, N_Time_0);

  assert(N == N_Time_0);

}
