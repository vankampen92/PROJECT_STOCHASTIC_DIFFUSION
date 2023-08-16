/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                             David Alonso, 2010 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

extern gsl_rng * r;   /* Global generator (define at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

// #define ASSERTION_TRUE

void Execute_One_Step(Community ** SP,
		      Parameter_Table * Table,
		      double max_Probability,
		      int * Event, int * x_Patch)
{
  int x, n_Event, n, n_Event_within_Age_Class, Q, k;
  Community * Patch;
  Parameter_Model * P = Table->P;

  Q = Table->TOTAL_No_of_DISEASE_STAGES * Table->TOTAL_No_of_AGE_CLASSES; /* Ex: 11 times 4 */

   /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  /* Hierarchic procedure to find the even to occur... */
  /* The event occurs in one of the local populations  */
  if(P->No_of_CELLS == 1) x = 0;
  else                                x = Choose_Village(max_Probability, SP, P);

  Patch = SP[x];  /* x represents the chosen patch undegoing a change. */

  n_Event = Discrete_Sampling(Patch->rToI, Table->TOTAL_No_of_EVENTS) - 1; /* 0, ..., 99 */

  n_Event_within_Age_Class = n_Event%Table->No_of_EVENTS;    /* 0, ..., 24 */
  n                        = n_Event/Table->No_of_EVENTS;    /* 0, ..., 3  */

  nS =  n0S  + n*Table->TOTAL_No_of_DISEASE_STAGES;   jS = nS + x*Q;
  nE =  n0E  + n*Table->TOTAL_No_of_DISEASE_STAGES;   jE = nE + x*Q;
  nI1=  n0I1 + n*Table->TOTAL_No_of_DISEASE_STAGES;  jI1 = nI1 + x*Q;
  nI2=  n0I2 + n*Table->TOTAL_No_of_DISEASE_STAGES;  jI2 = nI2 + x*Q;
  nA =  n0A  + n*Table->TOTAL_No_of_DISEASE_STAGES;   jA = nA + x*Q;
  nAd=  n0Ad + n*Table->TOTAL_No_of_DISEASE_STAGES;  jAd = nAd + x*Q;
  nY =  n0Y  + n*Table->TOTAL_No_of_DISEASE_STAGES;   jY = nY + x*Q;
  nR =  n0R  + n*Table->TOTAL_No_of_DISEASE_STAGES;   jR = nR + x*Q;
  aI =  a0I  + n*Table->TOTAL_No_of_DISEASE_STAGES;  jaI = aI + x*Q;
  aR =  a0R  + n*Table->TOTAL_No_of_DISEASE_STAGES;  jaR = aR + x*Q;
  aD =  a0D  + n*Table->TOTAL_No_of_DISEASE_STAGES;  jaD = aD + x*Q;

  switch( n_Event_within_Age_Class )
    {
    case  0:  /* (Age n): Infection (S --> S-1 and E --> E + 1) */                               /* 1 */
      Positivity_Control( 0, Table, x, jS, Y[jS], J[jS] ) ;
      Y[jS]--; J[jS]--;  Patch->n[nS]--;
      Y[jE]++; J[jE]++;  Patch->n[nE]++;

      break;
    case  1:  /* (Age n): Out-Migration (S --> S-1) and some other patch gains one */           /* Out S */
      Positivity_Control( 1, Table, x, jS, Y[jS], J[jS] );
      Y[jS]--; J[jS]--;  Patch->n[nS]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nS, Table);
      
      break;
    case  2:  /* 2 (Age n): In-Migration (S --> S+1) and some other patch loses one */         /* In S */
      Y[jS]++; J[jS]++;  Patch->n[nS]++;
      Local_Population_Increase(n,  Patch);
      Some_Other_Patch_Population_Decrease(x, n, nS, Table);
      
      break;
      
      
    case  3:  /* 3 (Age n):  Exposed into Infectious (E --> E-1 and I1 --> I1 + 1)*/             /* 2 */
      Positivity_Control( 3, Table, x, jE, Y[jE], J[jE]);
      Y[jE]--;  J[jE]--;   Patch->n[nE]--;
      Y[jI1]++; J[jI1]++;  Patch->n[nI1]++;
      

      break;
    case  4:  /* 4 (Age n): Out-Migration (E --> E-1) and some other patch gains one */        /* Out E */
      Positivity_Control( 4, Table, x, jE, Y[jE], J[jE]);
      Y[jE]--; J[jE]--;  Patch->n[nE]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nE, Table);

      break;
    case  5:  /* 5 (Age n): In-Migration (E --> E+1) and some other patch loses one */        /* In E */
       Y[jE]++; J[jE]++;  Patch->n[nE]++;
       Local_Population_Increase(n,  Patch);
       Some_Other_Patch_Population_Decrease(x, n, nE, Table);
      break;


    case  6:  /* 6 (Age n): Pre-Symptomatic into Infectious (I1 --> I1-1 and I2 --> I2 + 1) */   /* 3 */
              /*                                                         and AI --> AI + 1  */
      Positivity_Control( 6, Table, x, jI1, Y[jI1], J[jI1]);
      Y[jI1]--; J[jI1]--;  Patch->n[nI1]--;
      Y[jI2]++; J[jI2]++;  Patch->n[nI2]++;
      Y[jaI]++; J[jaI]++;  Patch->n[aI]++;

      break;
    case  7:  /* 7 (Age n): Pre-Symptomatic into A-symptomatic (I1 --> I1-1 and A --> A+1) */    /* 4 */
      Positivity_Control( 7, Table, x, jI1, Y[jI1], J[jI1]);
      Y[jI1]--; J[jI1]--;  Patch->n[nI1]--;
      Y[jA]++;  J[jA]++;   Patch->n[nA]++;

      break;
    case  8:  /* 8 (Age n): Out-Migration (I1 --> I1-1) and some other patch gains one */    /* Out I1  */
      Positivity_Control( 8, Table, x, jI1, Y[jI1], J[jI1]);
      Y[jI1]--; J[jI1]--;  Patch->n[nI1]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nI1, Table);

      break;
    case  9:  /* 9 (Age n): In-Migration (I1 --> I1+1) and some other patch loses one */     /* In I1 */
      Y[jI1]++; J[jI1]++;  Patch->n[nI1]++;
      Local_Population_Increase(n,  Patch);
      Some_Other_Patch_Population_Decrease(x, n, nI1, Table);

      break;


    case  10: /* 10 (Age n): Infectious I2 into Recovered (I2 -> I2-1 and R --> R + 1) */        /* 5 */
              /*                                                    and AR --> AR + 1  */
      Positivity_Control( 10, Table, x, jI2, Y[jI2], J[jI2]);
      Y[jI2]--; J[jI2]--;  Patch->n[nI2]--;
      Y[jR]++;  J[jR]++;   Patch->n[nR]++;
      Y[jaR]++; J[jaR]++;  Patch->n[aR]++;

      break;
    case  11: /* 11 (Age n): Infectious I2 into severe symptoms ( I2 -> I2-1 and Y --> Y + 1) */ /* 6 */
      Positivity_Control( 11, Table, x, jI2, Y[jI2], J[jI2]);
      Y[jI2]--; J[jI2]--;  Patch->n[nI2]--;
      Y[jY]++;  J[jY]++;   Patch->n[nY]++;

      break;
    case  12: /* 12 (Age n): Out-Migration (I2 --> I2-1) and some other patch gains one */      /* Out I2 */
      Positivity_Control( 12, Table, x, jI2, Y[jI2], J[jI2]);
      Y[jI2]--; J[jI2]--;  Patch->n[nI2]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nI2, Table);

      break;
    case  13: /* 13 (Age n): In-Migration (I2 --> I2 + 1) and some other patch loses one */    /* In I2 */
      Y[jI2]++; J[jI2]++;  Patch->n[nI2]++;
      Local_Population_Increase(n,  Patch);
      Some_Other_Patch_Population_Decrease(x, n, nI2, Table);

      break;


    case  14:  /* 14 (Age n): Assymptomatic into Recovered (A -> A-1 and R --> R + 1)    */
               /* We don't consider accumulating these recovered (AR --> AR + 1)
		  because they have not been detected by the system                      */      /* 7 */
      Positivity_Control( 14, Table, x, jA, Y[jA], J[jA]);
      Y[jA]--; J[jA]--;  Patch->n[nA]--;
			  Y[jR]++; J[jR]++;  Patch->n[nR]++;

      break;
    case  15:  /* 15 (Age n): Assymptomatic are detected   (A -> A-1 and Ad --> Ad + 1)  */      /* 8 */
               /*                                                    and AI --> AI + 1   */
      Positivity_Control( 15, Table, x, jA, Y[jA], J[jA]);
      Y[jA]--;  J[jA]--;   Patch->n[nA]--;
      Y[jAd]++; J[jAd]++;  Patch->n[nAd]++;
      Y[jaI]++; J[jaI]++;  Patch->n[aI]++;

      break;
    case  16:  /* 16 (Age n): Out-Migration (A --> A-1) and some other patch gains one */      /* Out A */
      Positivity_Control( 16, Table, x, jA, Y[jA], J[jA]);
      Y[jA]--; J[jA]--;  Patch->n[nA]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nA, Table);

      break;
    case  17:  /* 17 (Age n): In-Migration (A --> A +1) and some other patch loses one */     /* In A */
      Y[jA]++; J[jA]++;  Patch->n[nA]++;
      Local_Population_Increase(n,  Patch);
      Some_Other_Patch_Population_Decrease(x, n, nA, Table);

      break;


    case 18:  /* 18 (Age n): Detected Assymtomatic are recovered (Ad -> Ad-1 and R --> R + 1) */ /* 9 */
              /*                                                            and AR --> AR + 1 */
      Positivity_Control( 18, Table, x, jAd, Y[jAd], J[jAd]);
      Y[jAd]--; J[jAd]--;  Patch->n[nAd]--;
      Y[jR]++;  J[jR]++;   Patch->n[nR]++;
      Y[jaR]++; J[jaR]++;  Patch->n[aR]++;

      break;
    case  19:  /* 19 (Age n): Out-Migration (Ad --> Ad-1) and some other patch gains one */   /* Out Ad */
      Positivity_Control( 19, Table, x, jAd, Y[jAd], J[jAd] );
      Y[jAd]--; J[jAd]--;  Patch->n[nAd]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nAd, Table);

      break;
    case  20:  /* 20 (Age n): In-Migration (Ad --> Ad +1) and some other patch loses one */  /* In Ad */
      Y[jAd]++; J[jAd]++;  Patch->n[nAd]++;
      Local_Population_Increase(n,  Patch);
      Some_Other_Patch_Population_Decrease(x, n, nAd, Table);

      break;


    case  21:  /* 21 (Age n): Severe Infecious recover (Y --> Y-1 and R --> R + 1)    */        /* 10 */
               /*                                                 and AR --> AR + 1   */
               /*           : Severe Infectious don't move!!!                         */
      Positivity_Control( 21, Table, x, jY, Y[jY], J[jY] );
      Y[jY]--;  J[jY]--;   Patch->n[nY]--;
      Y[jR]++;  J[jR]++;   Patch->n[nR]++;
      Y[jaR]++; J[jaR]++;  Patch->n[aR]++;

      break;
    case  22:  /* 22 (Age n): Severe Infecious die (Y --> Y-1 and D --> D + 1 */                /* 11 */
      Positivity_Control( 22, Table, x, jY, Y[jY], J[jY] );
      Y[jY]--;  J[jY]--;   Patch->n[nY]--;
      Y[jaD]++; J[jaD]++;  Patch->n[aD]++;

      break;


    case  23:  /* 23 (Age n): Recovered move (R --> R-1 and some other patch gains one) */   /* Out R */
      Positivity_Control( 23, Table, x, jR, Y[jR], J[jR]);
      Y[jR]--; J[jR]--;  Patch->n[nR]--;
      Local_Population_Decrease(n,  Patch);
      Some_Other_Patch_Population_Increase(x, n, nR, Table);

      break;
    case  24:  /* 24 (Age n): In-Migration (R --> R + 1) and some other patch loses one */   /* In R */
      Y[jR]++; J[jR]++;  Patch->n[nR]++;
      Local_Population_Increase(n,  Patch);
      Some_Other_Patch_Population_Decrease(x, n, nR, Table);;

      break;

    default:
    /* Something is very very wrong!!! */
      printf("The number of event occurring should be between 0 and 24\n");
      printf("Event to Occur = %d\n", n_Event);
      Print_Press_Key(1,0,".");
      exit(0);
    }


  (*Event) = n_Event;  (*x_Patch) = x;
}

void Local_Population_Decrease ( int n, Community * Patch )
{

  switch( n )
    {
    case  0:  Patch->N0--;
      break;
    case  1:  Patch->N1--;
      break;
    case  2:  Patch->N2--;
      break;
    case  3:  Patch->N3--;
      break;
    default:
      printf("Age class out of range\n");
      printf("Age class can be 0, 1, 2, or 4, but Age Class = %d\n", n);
      printf("The program will exit\n");
      exit(0);
    }
}

void Local_Population_Increase ( int n, Community * Patch )
{

  switch( n )
    {
    case  0:  Patch->N0++;
      break;
    case  1:  Patch->N1++;
      break;
    case  2:  Patch->N2++;
      break;
    case  3:  Patch->N3++;
      break;
    default:
      printf("Age class out of range\n");
      printf("Age class can be 0, 1, 2, or 4, but Age Class = %d\n", n);
      printf("The program will exit\n");
      exit(0);
    }
}

void Some_Other_Patch_Population_Decrease(int x, int a, int nS,
					  Parameter_Table * Table)
{
  /* Input:

     .  x: Patch label
     .  a: Age Class;
     . nS: Diseasea Status  (0, ..., 43)
     . Table: Parameter Table

  */
  int Q, i, j, n, n_Patch;
  Community ** Patch = Table->Patch_System;

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  Q = Table->TOTAL_No_of_DISEASE_STAGES * Table->TOTAL_No_of_AGE_CLASSES; /* Ex: 11 times 4 */

  double * Imm_Preassure_Vector = (double *)calloc( Patch[x]->No_NEI, sizeof(double) );

  for(n=0; n<Patch[x]->No_NEI; n++)
    Imm_Preassure_Vector[n] = Patch[x]->In_Migration_Vector[a][n]*(double)Patch[x]->NEI[n]->n[nS];
    
  n_Patch = Discrete_Sampling( Imm_Preassure_Vector, Patch[x]->No_NEI ) - 1;
  
  assert(Patch[Patch[x]->Patch_Connections[n_Patch]] == Patch[x]->NEI[n_Patch]);
  
  j = nS + Patch[x]->Patch_Connections[n_Patch]*Q;
  
  Positivity_Control( 100, Table, Patch[x]->Patch_Connections[n_Patch], j, Y[j], J[j] );
  Y[j]--; J[j]--;  Patch[x]->NEI[n_Patch]->n[nS]--;
  Local_Population_Decrease ( a, Patch[x]->NEI[n_Patch] );

  free(Imm_Preassure_Vector); 
}

void Some_Other_Patch_Population_Increase(int x, int a, int nS,
					  Parameter_Table * Table)
{
  /* Input:

     .  x: Patch label
     .  a: Age Class;
     . nS: Diseasea Status  (0, ..., 43)
     . Table: Parameter Table

  */
  int Q, k, j, n_Patch;
  Community ** Patch = Table->Patch_System;

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  Q = Table->TOTAL_No_of_DISEASE_STAGES * Table->TOTAL_No_of_AGE_CLASSES; /* Ex: 11 times 4 */

  n_Patch = Discrete_Sampling(Patch[x]->Out_Migration_Vector[a], Patch[x]->No_NEI) - 1;

  j = nS + Patch[x]->Patch_Connections[n_Patch]*Q;

  Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[nS]++;
  Local_Population_Increase ( a, Patch[x]->NEI[n_Patch] );

}

void Positivity_Control( int Event, Parameter_Table * Table,
			 int x, int jS, double Y, int J)
{
  int Q, nS, Non_Positive;
  Community ** Patch = Table->Patch_System;

#if defined ASSERTION_TRUE

  Q = Table->TOTAL_No_of_DISEASE_STAGES * Table->TOTAL_No_of_AGE_CLASSES; /* Ex: 11 times 4 */
  nS = jS%Q;

  Non_Positive = 0;

  if ( Y <= 0.0 )            Non_Positive = 1;

  if ( J <= 0 )              Non_Positive = 1;

  if ( Patch[x]->n[nS] <= 0) Non_Positive = 1;

  if( Non_Positive == 1 ) {
    printf (" Event No %d:", Event);  
    printf (" Y[%s] = %g\t", Table->Model_Variable_Symbol[jS], Y);
    printf ("J[%s] = %d\t",  Table->Model_Variable_Symbol[jS], J);
    printf ("n[%s] = %d\n",  Table->Model_Variable_Symbol[jS], Patch[x]->n[nS]);
    exit(0);
  }
  
#endif
}
