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
  if(P->No_of_CELLS == 1) x = y = 0;
  else                    x = y = Choose_Village(max_Probability, SP, P);

  /* x and y will end up differing only in case there is movemnt event!!! */

  Patch = SP[x];  /* x represents the chosen patch undegoing a change. */

  if(Table->TOTAL_No_of_EVENTS > 1)
    n_Event = Discrete_Sampling(Patch->rToI, Table->TOTAL_No_of_EVENTS) - 1; /* 0, ..., 10 */
    
  else {
    printf(" The total number of events tant potentially could happen patch %d\n", x);
    printf(" is zero??? (TOTAL_No_of_EVENTS = %d)\n", Table->TOTAL_No_of_EVENTS); 
    printf(" Something very wrong with your code\n");
    printf(" The program will exit\n");
    Print_Press_Key(1,0,"."); 
    exit(0); 
  }

  R   = x*Table->LOCAL_STATE_VARIABLES + Table->R;
  A   = x*Table->LOCAL_STATE_VARIABLES + Table->A;
  RA  = x*Table->LOCAL_STATE_VARIABLES + Table->RA;
  ARA = x*Table->LOCAL_STATE_VARIABLES + Table->ARA;

  assert( n_Event < Table->TOTAL_No_of_EVENTS );
  
  // Print_Meta_Community_Patch_System (Table);
  
  switch( n_Event )
    {
    case  0:  /* Out-Migration (R --> R-1) from patch x and some other patch gains one *//* Out R */
      Positivity_Control( 0, Table, x, R, Y[R], J[R] );

      Y[R]--; J[R]--; Patch->n[Table->R]--;

      y = Some_Other_Patch_Population_Increase(x, Table->R, Table);

      break;
      
    case  1:  /* External Immigration from outside the system (R --> R+1) *//* External Imm R */

      Y[R]++; J[R]++;  Patch->n[Table->R]++;
      
      break;
      
    case  2:  /* Local Death  (R --> R-1)  */                               /* Death R */
      Positivity_Control( 2, Table, x, R, Y[R], J[R] );

      Y[R]--; J[R]--;  Patch->n[Table->R]--;

      break; 

    case  3:  /* Out-Migration (A --> A-1) from patch x and some other patch gains one *//* Out A */
      Positivity_Control( 3, Table, x, A, Y[A], J[A] );

      Y[A]--; J[A]--; Patch->n[Table->A]--;

      y = Some_Other_Patch_Population_Increase(x, Table->A, Table);

      break;
      
    case  4:  /* External Immigration from outside the system (A --> A+1) *//* External Imm A */

      Y[A]++; J[A]++;  Patch->n[Table->A]++;
      
      break;
      
    case  5:  /* Local Death  (A --> A-1)  */                               /* Death A */
      Positivity_Control( 5, Table, x, A, Y[A], J[A] );

      Y[A]--; J[A]--;  Patch->n[Table->A]--;

      break; 

    case 6: /* Consumer Consumption of Resource and dimmer formation */ /* R + A ---> RA */
      Positivity_Control( 6, Table, x, R, Y[R], J[R] );
      Positivity_Control( 6, Table, x, A, Y[A], J[A] );

      Y[R]--; J[R]--; Patch->n[Table->R]--;
      Y[A]--; J[A]--;  Patch->n[Table->A]--;
      
      Y[RA]++; J[RA]++;  Patch->n[Table->RA]++;
      break;

    case 7: /* Dimmer degratation into two individual consumers (Consumer Growth) */ /* RA ---> A + A */
      Positivity_Control( 7, Table, x, RA, Y[RA], J[RA] );
      
      Y[A]++; J[A]++;  Patch->n[Table->A]++;
      
      Y[A]++; J[A]++;  Patch->n[Table->A]++;
      
      Y[RA]--; J[RA]--;  Patch->n[Table->RA]--;
      
      break;

    case 8: /* Local Growth of Resources  */                       /* R  ---> R + R */

      Y[R]++; J[R]++;  Patch->n[Table->R]++;
      
      break;

    case 9: /* Consumer Interference */                           /* RA + A ---> ARA */
      Positivity_Control( 9, Table, x, RA, Y[RA], J[RA] );
      Positivity_Control( 9, Table, x, A, Y[A], J[A] );
      
      Y[RA]--; J[RA]--; Patch->n[Table->RA]--;
      Y[A]--; J[A]--;  Patch->n[Table->A]--;
      
      Y[ARA]++; J[ARA]++;  Patch->n[Table->ARA]++;
      
      break;

    case 10: /* Degradation of triplets */                        /* ARA ---> RA + A */
      Positivity_Control( 10, Table, x, ARA, Y[ARA], J[ARA] );

      Y[RA]++; J[RA]++; Patch->n[Table->RA]++;
      Y[A]++; J[A]++;  Patch->n[Table->A]++;
      
      Y[ARA]--; J[ARA]--;  Patch->n[Table->ARA]--;
      
      break;
        
    default:
    /* Something is very very wrong!!! */
      printf("The number of event occurring should be between 0 and 0\n");
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

  j = Sp + Patch[x]->Patch_Connections[n_Patch]*Table->LOCAL_STATE_VARIABLES;

  Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;

  return(Patch[x]->Patch_Connections[n_Patch]);
}

void Positivity_Control( int Event, Parameter_Table * Table,
			 int x, int jS, double Y, int J)
{
  int Q, nS, Non_Positive;
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
    exit(0);
  }

#endif
}
