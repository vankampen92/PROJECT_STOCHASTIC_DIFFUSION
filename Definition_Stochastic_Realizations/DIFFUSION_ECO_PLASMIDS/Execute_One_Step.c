/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                             David Alonso, 2010 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

extern gsl_rng * r;   /* Global generator variable defined at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

void Calculate_DonorID_and_RecipientID(Parameter_Table * Table, int n, 
                                       int * Donor_ID, int * Recipient_ID);
int Calculate_Transconjugant_Strain_ID_from_Profile(Parameter_Table * Table, 
                                                    int Sp_Strain, int * Profile );                                       
int Conjugation_Transconjugant_Calculation(Parameter_Table * Table, 
                                           int x, int Sp, int Sp_1, 
                                           int * Sp_2);

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
  /* Calling hierarchy: 
     main() (main.c) 
      --->  M_O_D_E_L___S_T_O( &Table ); (MODEL_STO.c) 
        --->  S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S (...); (./Definition_Stochastic_Realizations/Stochastic_Time_Dynamics.c)
          ---> Advance_Current_Time(...) (./Definition_Stochastic_Realizations/Advance_Current_Time.c)
            ---> Exectute_One_Step(...) (./Definition_Stochastic_Realizations/DIFFUSION_ECO_PLASMIDS/Execute_One_Step.c) 
  */

  int x, y, n_Event, Index_Patch, Index_Event, Index_Local, n, n_Event_Sp, Sp, j;
  Community * Patch;
  Parameter_Model * P = Table->P;

  int Sp_1, Sp_2, R_1, R_2, Sp_Strain, Sp_Profile, boolean_Success;

  /* These are bacterial strains IDs. They are intialized with an impossible FLAG */
  /* because 0 <= Sp_1 < Table->No_of_RESOURCES-1                                */
  Sp_1 = Table->No_of_RESOURCES; 
  Sp_2 = Table->No_of_RESOURCES;
  boolean_Success = 1; /* Always One Step Exectud (Success)!!! If failure, this variable is turned into 0 */ 

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
     x and y will end up being different only in case there is a movemnt event!!! 
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
  
  // Sp is the FOCAL SPECIES undergoing the event n_Event_Sp
  if(n_Event < (Table->No_of_EVENTS-1)*Table->No_of_RESOURCES){ 
    Sp       = n_Event/(Table->No_of_EVENTS-1); /* Sp should be an ID from 0 to Table->No_of_RESOURCES -1 */
    n_Event_Sp = n_Event%(Table->No_of_EVENTS-1); /* Event to happen within Sp species                    */
  }
  else{
    n_Event_Conjugation = n_Event-(Table->No_of_EVENTS-1)*Table->No_of_RESOURCES;
    Calculate_DonorID_and_RecipientID(Table, n_Event_Conjugation, &Sp, &Sp_1);
    /* More efficient just looking at a table... rather than calculating them each time !!! */
    int SSp   = Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[n_Event_Conjugation][0]; /* DONOR     */
    int SSp_1 = Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[n_Event_Conjugation][1]; /* RECIPIENT */

    assert(SSp == Sp && SSp_1 == Sp_1);
    n_Event_Sp = 6; 
  }

  R   = x*Table->LOCAL_STATE_VARIABLES + Sp;      

  assert( Sp >= 0   && Sp < Table->No_of_RESOURCES );
  assert( Sp_1 >= 0 && Sp_1 < Table->No_of_RESOURCES );
  assert( n_Event < Table->TOTAL_No_of_EVENTS );
  assert( x       < Table->No_of_CELLS );

  // Print_Meta_Community_Patch_System (Table);

  switch( n_Event_Sp )
    {
    case  0:  /* 0: Bacteria Out-Migration (R --> R-1) and some other patch gains one */ /* Outmigraiton (R) */
      Positivity_Control( 0, Table, x, Sp, Y[R], J[R] ); 
      Y[R]--; J[R]--; Patch->n[Sp]--;                          Patch->m_0++; 
                      Patch->Local_Strain_Population[Sp]->n--;

      y = Some_Other_Patch_Population_Increase(x, Sp, Table); /* Sp increases in patch y where alos empty space 
                                                                 decreases */
      break;

    case  1:  /* External Immigration on focal patch x from outside the system (R --> R+1) */ /* External Imm R */
      Carrying_Capacity_Control( 1, Table, x, Sp, Y[R], J[R] );

      Y[R]++; J[R]++;  Patch->n[Sp]++;                         Patch->m_0--; 
                       Patch->Local_Strain_Population[Sp]->n++;       
    break;

    case  2:  /* Local Death of Propagules (R --> R-1)  */                                /* Death (R) */
      Positivity_Control( 2, Table, x, Sp, Y[R], J[R] );
      Y[R]--; J[R]--;  Patch->n[Sp]--;                          Patch->m_0++; 
                       Patch->Local_Strain_Population[Sp]->n--;
      break;

    case  3:  /* 3: Competition induced mortality on Species R  (R --> R-1) */            /* Production (R) */
      Positivity_Control( 3, Table, x, Sp, Y[R], J[R] );
      Y[R]--; J[R]--;  Patch->n[Sp]--;                          Patch->m_0++; 
                       Patch->Local_Strain_Population[Sp]->n--;
      break;

    case  4:  /* 4: Cell dividion (without error) */                       /* Cell Division withot error (R) */
      Carrying_Capacity_Control( 4, Table, x, Sp, Y[R], J[R] );
      Y[R]++;  J[R]++;  Patch->n[Sp]++;                         Patch->m_0--;   /* Free space gets reduced */
                        Patch->Local_Strain_Population[Sp]->n++;
      break;

    case 5:   /* 5: Cell division (with error) */                        /* Cell Division with error */   
      Calculate_Strain_and_Profile(Table, Sp, &Sp_Strain, &Sp_Profile);
      
      int Sp_Strain_0  = Table->StrainType_and_Profile[Sp][0];     /* Strain Type          */
      int Sp_Profile_0 = Table->StrainType_and_Profile[Sp][1];     /* Profile No Specifier */

      assert(Sp_Strain_0 == Sp_Strain && Sp_Profile_0 == Sp_Profile);
      
      Sp_1 = Table->n_0[Sp_Strain];                                      /* Strain ID of the plasmid free type */
      R_1  = x*Table->LOCAL_STATE_VARIABLES + Sp_1;

      Carrying_Capacity_Control( 5, Table, x, Sp_1, Y[R_1], J[R_1] );
      Y[R_1]++;  J[R_1]++;  Patch->n[Sp_1]++;                   Patch->m_0--;   /* Free space gets reduced */
                            Patch->Local_Strain_Population[Sp_1]->n++;
      break; 

    case 6:  /* 6: Cell to Cell Conjugation  */  
             /* Donor(k, p) + Recipient(l, q) --->  Donor(k, p) + Transconjugant(l, r) */
             /* Three species IDs may be involved: 
                Sp, Donor(k, p); 

                R, R_1, R_2 should stay in the range:
                [0, Table->No_of_CELLS*Table->LOCAL_STATE_VARIABLES-1]
                
                Sp_1 = Recipient(l, q); 
                Sp_2 = Transconjugant(l, r) 
             */

      boolean_Success = Conjugation_Transconjugant_Calculation(Table, x, Sp, Sp_1, &Sp_2); 
      if (boolean_Success == 1) { /* Effective Conjugation */
        /* Recipient decreases by 1 */
        R_1 = x*Table->LOCAL_STATE_VARIABLES + Sp_1;

        Positivity_Control( 6, Table, x, Sp_1, Y[R_1], J[R_1] );        
        Y[R_1]--;  J[R_1]--;  Patch->n[Sp_1]--;                       /* Free space remains the same */
                              Patch->Local_Strain_Population[Sp_1]->n--;
        
        /* Transconjugant increases by 1 */
        R_2 = x*Table->LOCAL_STATE_VARIABLES + Sp_2;
        
        Carrying_Capacity_Control( 6, Table, x, Sp_2, Y[R_2], J[R_2] );
        Y[R_2]++;  J[R_2]++;  Patch->n[Sp_2]++;                       /* Free space remains the same */
                              Patch->Local_Strain_Population[Sp_2]->n++;
      }
      else {                      /* Conjugation has failed */
        /* No configurational change in the system: the system remains the same and
           no update of the rates is required 
        */
        Sp_1 = Sp;  /* This is a conjugation failure */
        Sp_2 = Sp;
        boolean_Success = 0;   
      }
      break;
      
    default:
    /* Something is very very wrong!!! */
      printf("The number of event occurring should be between 0 and 12\n");
      printf("Event to Occur = %d\n", n_Event);
      Print_Press_Key(1,0,".");
      exit(0);
    }

  (*Event) = n_Event;  /* From n_Event, the focal species can be decoded */  
  x_Patch[0] = x;  
  x_Patch[1] = y;  
  x_Patch[2] = Sp_1;   /* RECIPIENT */
  x_Patch[3] = Sp_2;   /* DONOR     */
  x_Patch[4] = boolean_Success; 
  
  /* (*Event) allows us to identify both the affected species ID 
     (donor, if the event has been a conjugation) and the patch where 
     the event has occured, x */
  /* If the event is a between-patch movement, then x is the affected patch, 
     and y is the patch where the individual bacterial cell has moved to from patch x. 
     Otherwise, always x == y. 
  */
  /* If conjugation has occurred, Sp_1 is the ID of recipient species, 
     and Sp_2, the ID of the transconjugant species */ 
}

int Some_Other_Patch_Population_Increase(int x, int Sp,
					                               Parameter_Table * Table)
{
  /* Input:

     . x: Patch label of the patch from which the individuals goes out
     . Sp: Species;
     . Table: Parameter Table

     Output:

     . y is Patch[x]->Patch_Connections[n_Patch], i.e., the Patch label 
       of the patch receiving the immigrant

  */
  int k, j, y, n_Patch;
  Community ** Patch = Table->Patch_System;
  int R; 

  int    * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  n_Patch = Discrete_Sampling(Patch[x]->Out_Migration_Vector[Sp], Patch[x]->No_NEI) - 1;

  j = Sp + Patch[x]->Patch_Connections[n_Patch] * Table->LOCAL_STATE_VARIABLES;

  y = Patch[x]->Patch_Connections[n_Patch];
  R = y * Table->LOCAL_STATE_VARIABLES + Sp;
  Carrying_Capacity_Control( 0, Table, y, Sp, Y[R], J[R] );
  
  Y[j]++; J[j]++;  Patch[x]->NEI[n_Patch]->n[Sp]++;
                   Patch[x]->NEI[n_Patch]->Local_Strain_Population[Sp]->n++;

  return(Patch[x]->Patch_Connections[n_Patch]);
}

void Calculate_DonorID_and_RecipientID(Parameter_Table * Table, int n, 
                                       int * Donor_ID, int * Recipient_ID)
{
  /* Input: 
     . n, label of the conjugation event
     Output:
     . Donor_ID, from 0 to No_of_RESOURCES-1.
     . Recipient_ID, from 0 to No_of_RESOURCES-1.
  */
  int i, j;
  int N; 

  N = Table->n_R[0]; 
  i = 0; 
  while ( (n - N) >= 0 )
  {
    N += Table->n_R[i+1];
    i++;   
  }
  * Donor_ID = i; 
  
  N = 0; 
  for(i=0; i < (* Donor_ID); i++)
  { 
    N += Table->n_R[i]; 
  }
  j = n - N;
  * Recipient_ID = Table->Recipient_List_Indeces[* Donor_ID][j]; 
}

int Conjugation_Transconjugant_Calculation(Parameter_Table * Table, 
                                           int x, int Sp, int Sp_1, 
                                           int * Sp_2)
{
  int Success; 
  double S; 
  int i, j, k, n, m, m_0;
  int Sp_Strain, Sp_Profile; 
  Community ** Patch = Table->Patch_System;
  /* Sp:       Donor ID          */
  /* Sp_1:     Recipient ID      */
  /* (* Sp_2): Transconjugant ID */
  
  /* Calculation of the Transconjugant  */
  int * Transconjugant_Profile = (int *)calloc(Table->No_of_PLASMIDS, 
                                               sizeof(int));
  /* ------------------------------------- */  
  /* 1: Counting Plasmids in the Donor and 
        Determining the identity and number of plasmids that 
        will be transfered (n) 
  */
  int    * Pl = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int));
  n = 0; S = 1.0; 
  for(k = 0; k<Table->No_of_PLASMIDS; k++) 
    if (PATCH[x]->Local_Strain_Population[Sp]->Profile[k] == 1) { 
      if( (RANDOM) < Table->Chi_C_0 ) {
        Pl[n] = k;
        n++;
      } 
    } 
  
  /* The conjugation event will be avorted according to the philosophy 
     or "all or nothing", this is, in the case either at least one plasmid 
     in the potential transfering set cannot infect the recipient cell or at 
     least there is one single incompatibility between any of the potentially 
     transfered plasmids and the plasmids present in the recipient cell 
     (then No Success, i.e., Success = 0). 
  */
  m = 0; /* Incompatibility Counter */
  /* 2: Assessing if this set of 'n' plasmids can all infect the recipient */
  Calculate_Strain_and_Profile(Table, Sp_1, &Sp_Strain, &Sp_Profile);
  
  int Sp_Strain_0  = Table->StrainType_and_Profile[Sp_1][0];     /* Strain Type          */
  int Sp_Profile_0 = Table->StrainType_and_Profile[Sp_1][1];     /* Profile No Specifier */

  assert(Sp_Strain_0 == Sp_Strain && Sp_Profile_0 == Sp_Profile);

  for(i=0; i<n; i++)  
    if(Table->IBP[Sp_Strain][Pl[i]] == 0.0)
      m++;

  if( m == 0) {
  /* 3: Assessing that at least one of the transfered plasmids was not already
        present in the recipient and building the potential tranconjugant 
        profile 
  */
    m_0 = 0; 
    for(i = 0; i<n; i++) 
      if (PATCH[x]->Local_Strain_Population[Sp_1]->Profile[Pl[i]] == 1) 
        m_0++;

    if( m_0 == n) /* All transfered plasmids were already present */
      m++; 
  }

  if( m == 0) {
    /* 4: Assessing plasmid-plasmid incompatibilities in the transconjugant */
    /* Ix[] will store the list of plasmids in the recipient, while 
       Pl[] already stores the list of plasmids to be transfered 
    */ 
    int    * Ix = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int)); 
    m_0 = 0;
    for(i = 0; i < Table->No_of_PLASMIDS; i++) 
      if (Patch[x]->Local_Strain_Population[Sp_1]->Profile[i] == 1)
        Ix[m_0++] = i;  

    /* m_0 plasminds were already present in the recipient */
    /* n plasmids will be transfered */
    for(i = 0; i < m_0; i++) {
      Tranconjugant_Profile[Ix[i]] = 1; 
      for(j = 0; j < n; j++){
        if (Table->CPP[Ix[i]][Pl[j]] == 0.0) m++;

        Tranconjugant_Profile[Pl[j]] = 1; 
      }
    }

    free(Ix);
  }
    
  if( m == 0)
  { 
    /* Looking for the transconjugant ID in the different profiles associated 
       to recipient strain ID (* Sp_1) with strain type Sp_Strain 
    */
    * Sp_2 = Calculate_Transconjugant_Strain_ID_from_Profile(Table, Sp_Strain, 
                                                             Transconjugant_Profile); 

    if( (* Sp_2) == Table->No_of_RESOURCES ) {

      printf(" No transconjunt profile found!!\n");
      printf(" Some inconsistency in the construction of the transconjugant!!!");
      printf(" Error in: Conjugation_Transconjugant_Calculation()\n");

      exit(0);
    }
    
    Success = 1;
  } 
  else 
    Success = 0;   

  free(Pl); 
  free(Transconjugant_Profile);
  
  return (Success); /* Success = 1, conjugation success  
                       Success = 0, conjugation failure 
                    */
} 

int Calculate_Transconjugant_Strain_ID_from_Profile(Parameter_Table * Table, 
                                                    int Sp_Strain, int * Profile )
{
  /* Input:
     . Table
     . Sp_Strain, Strain Type: Sp_Strain >= 0 && Sp_Strain < No_of_STRAINS
     . Profile, Potential profile to search for. 

     Output:
     . Strain ID corresponding to input Strain Type 'Sp_Strain' and 'Profile'  
  */
  
}

