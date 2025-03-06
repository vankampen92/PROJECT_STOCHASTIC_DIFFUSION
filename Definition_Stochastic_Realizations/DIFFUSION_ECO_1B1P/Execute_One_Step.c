/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                             David Alonso, 2025 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

extern gsl_rng * r;   /* Global generator variable defined at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

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
            ---> Exectute_One_Step(...) (./Definition_Stochastic_Realizations/DIFFUSION_ECO_1B1P/Execute_One_Step.c) 
  */
  int Conjugation_Success; 
  int x, y, n_Event, Index_Patch, Index_Event, Index_Local, n, j, R_0, R_1;
  Community * Patch;
  Parameter_Model * P = Table->P;

  int Sp_0, Sp_1, boolean_Success;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  int * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  /* Hierarchic procedure to find the even to occur... */
  /* The event occurs in one of the local populations  */
  if (P->No_of_CELLS == 1) {
    x = y = 0;
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Choose_Village_and_Event_Binary_Tree(max_Probability, SP, P, 
                                           &Index_Patch, &Index_Event);
      x = y = Index_Patch; 
      n_Event = Index_Event;
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Choose_Village_and_Event_Priority_Queu(max_Probability, SP, P, 
                                             &Index_Patch, &Index_Event);
      x = y = Index_Patch; 
      n_Event = Index_Event;

      assert(Index_Event == Table->Treeroot->index);
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
      x = y = Index_Patch; 
      n_Event = Index_Event;
    #elif defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Choose_Village_and_Event_Priority_Queu(max_Probability, SP, P, 
                                             &Index_Patch, &Index_Event);
      x = y = Index_Patch; 
      n_Event = Index_Event;

      assert((Index_Patch * Table->TOTAL_No_of_EVENTS + Index_Event) == Table->Treeroot->index);
      assert(Table->Tree_Node_Index[Table->Treeroot->index] == Table->Treeroot);
    #else
      x = y = Choose_Village(max_Probability, SP, P);
    #endif
  }

  /* When this function 'returns', 
     x and y will end up being different only in case there is a movemnt event!!! 
  */
  Patch = SP[x];  /* x represents the chosen patch undergoing a change. */
  
  #if defined EVENT_DISCRETE_SAMPLING
    if (Table->TOTAL_No_of_EVENTS > 1) 
      n_Event = Discrete_Sampling(Patch->rToI, Table->TOTAL_No_of_EVENTS) - 1; 
      /* 0, ..., Table->TOTAL_No_of_EVENTS-1 */
  #endif
  
  if (Table->TOTAL_No_of_EVENTS == 0) {
    printf(" The potential number of events that potentially could happen in patch %d\n", x);
    printf(" is zero??? (TOTAL_No_of_EVENTS = %d)\n", Table->TOTAL_No_of_EVENTS);
    printf(" If it is, there is something very wrong with your code\n");
    printf(" The program will exit!!!\n");
    Print_Press_Key(1,1,"Error in Execute_One_Step() function");
    exit(0);
  }

  // Print_Discrete_Probability_Distribution(Table, n_Event, x); /* Comment out if it works */

  if (n_Event < 0 || n_Event >= Table->TOTAL_No_of_EVENTS) {
    printf("n_Event = %d\n", n_Event);
    printf("This event is not defined\n");
    printf("The program will exit\n");
    exit(0);
  }
 
  R_0 = x * Table->LOCAL_STATE_VARIABLES + 0; /* A_0 represents the plasmid-free cell (or population) */
  R_1 = x * Table->LOCAL_STATE_VARIABLES + 1; /* A_1 represents plasmid-carrying cell (or population) */   

  assert(x < Table->No_of_CELLS && x >= 0);

  // Print_Meta_Community_Patch_System (Table);
  Conjugation_Success = 0;
  
  switch (n_Event) {
    case 0:  /* 0: Bacteria Out-Migration (A_0 --> A_0 - 1) and some other patch gains one */ /* Outmigration (A_0) */
      assert(Table->Mu_RP[0] > 0.0 && Table->No_of_CELLS > 0);

      Positivity_Control(0, Table, x, 0, Y[R_0], J[R_0]); 
      Y[R_0]--; J[R_0]--; Patch->n[0]--; Patch->m_0++;         // Patch->Local_Strain_Population[0]->n--;

      y = Some_Other_Patch_Population_Increase(x, 0, Table); /* A_0 increases in patch y where also empty space 
                                                                decreases */

      if( y == x ) /* No movement can take place. Only the loss in the focal patch x */
        n_Event = 2; /* Local Death of Plasmid Free Cells (A_0 --> A_0 - 1)  */ /* Death (A_0) */   
                                                              
      x_Patch[2] = 0;  x_Patch[3] = 0;                                                             
      break;

    case 1:  /* 1: External Immigration on focal patch x from outside the system (A_0 --> A_0 + 1) */ /* External Imm A_0 */
      Carrying_Capacity_Control(1, Table, x, 0, Y[R_0], J[R_0]);
      Y[R_0]++; J[R_0]++; Patch->n[0]++; Patch->m_0--;          // Patch->Local_Strain_Population[0]->n++;       
      
      x_Patch[2] = 0;  x_Patch[3] = 0;      
      break;

    case 2:  /* 2: Local Death of Plasmid Free Cells (A_0 --> A_0 - 1)  */                       /* Death (A_0) */
      Positivity_Control(2, Table, x, 0, Y[R_0], J[R_0]);
      Y[R_0]--; J[R_0]--; Patch->n[0]--; Patch->m_0++;           // Patch->Local_Strain_Population[0]->n--;
      
      x_Patch[2] = 0;  x_Patch[3] = 0;
      break;

    case 3:  /* 3: Competition induced mortality on Species A_0  (A_0 --> A_0 - 1) */            /* Death (A_0) */
      Positivity_Control(3, Table, x, 0, Y[R_0], J[R_0]);
      Y[R_0]--; J[R_0]--; Patch->n[0]--; Patch->m_0++;           // Patch->Local_Strain_Population[0]->n--;

      x_Patch[2] = 0;  x_Patch[3] = 0;
      break;

    case 4:  /* 4: Cell division (without error) */                       /* Cell Division without error (A_0 --> A_0 + 1) */
      Carrying_Capacity_Control(4, Table, x, 0, Y[R_0], J[R_0]);
      Y[R_0]++; J[R_0]++; Patch->n[0]++; Patch->m_0--;           // Patch->Local_Strain_Population[0]->n++;

      x_Patch[2] = 0;  x_Patch[3] = 0;
      break;

    case 5:  /* 5: Bacteria Out-Migration (A_1  --> A_1  - 1) and some other patch gains one */ /* Outmigration (A_1) */   
      assert(Table->Mu_RP[1] > 0.0 && Table->No_of_CELLS > 0);

      Positivity_Control(5, Table, x, 1, Y[R_1], J[R_1]); 
      Y[R_1]--; J[R_1]--; Patch->n[1]--; Patch->m_0++;            // Patch->Local_Strain_Population[1]->n--;

      y = Some_Other_Patch_Population_Increase(x, 1, Table); /* Sp increases in patch y where 
                                                                also empty space decreases */
      
      if( y == x ) /* No movement can take place. Only the loss in the focal patch x */
        n_Event = 7; /* Local Death of Plasmid Carrying Cells (A_1 --> A_1 - 1)  */ /* Death (A_1) */   
      
      x_Patch[2] = 1;  x_Patch[3] = 1;                                                          
      break;

    case 6:  /* 6: External Immigration on focal patch x from outside the system (A_1  --> A_1  + 1) */ /* External Imm A_1 */
      Carrying_Capacity_Control(6, Table, x, 1, Y[R_1], J[R_1]);

      Y[R_1]++; J[R_1]++; Patch->n[1]++; Patch->m_0--;             // Patch->Local_Strain_Population[1]->n++; 
      
      x_Patch[2] = 1;  x_Patch[3] = 1;                                                          
      break;

    case 7:  /* 7: Local Death of Plasmid-Carrying Cells (A_1  --> A_1  - 1)  */                         /* Death (A_1) */
      Positivity_Control(7, Table, x, 1, Y[R_1], J[R_1]);
      Y[R_1]--; J[R_1]--; Patch->n[1]--; Patch->m_0++;             // Patch->Local_Strain_Population[1]->n--;
      
      x_Patch[2] = 1;  x_Patch[3] = 1;                                                          
      break;

    case 8:  /* 8: Competition induced mortality on Species R  (A_1  --> A_1 - 1) */                   /* Death (A_1) */
      Positivity_Control(8, Table, x, 1, Y[R_1], J[R_1]);
      Y[R_1]--; J[R_1]--; Patch->n[1]--; Patch->m_0++;             // Patch->Local_Strain_Population[1]->n--;
      
      x_Patch[2] = 1;  x_Patch[3] = 1;
      break;

    case 9:  /* 9: Cell division (without error)  A_1  ---> A_1 + 1 */                       /* Cell Division without error (A_1  --> A_1  + 1) */
      Carrying_Capacity_Control(9, Table, x, 1, Y[R_1], J[R_1]);
      Y[R_1]++; J[R_1]++; Patch->n[1]++; Patch->m_0--;             // Patch->Local_Strain_Population[1]->n++;
      
      x_Patch[2] = 1;  x_Patch[3] = 1;
      break;

    case 10:   /* 10: Cell division (with error)  A_0 ---> A_0 + 1 */                        /* Cell Division with error (A_1 ---> A_1  + 1) */   
      Carrying_Capacity_Control(10, Table, x, 0, Y[R_0], J[R_0]);
      Y[R_0]++; J[R_0]++; Patch->n[0]++; Patch->m_0--;             // Patch->Local_Strain_Population[0]->n++;
      
      x_Patch[2] = 0;  x_Patch[3] = 0;
      break;            
    
    case 11:  /* 11: Conjugation: A_1 + A_0 ---> A_1 + A_1 */
      /* Donor(k, p) + Recipient(l, q) --->  Donor(k, p) + Transconjugant(l, r) */
      /* The event occurs in the focal patch x: Transconjugant increases */
      Carrying_Capacity_Control(11, Table, x, 1, Y[R_1], J[R_1]);
      Y[R_1]++; J[R_1]++; Patch->n[1]++;                                                 
                         // Patch->Local_Strain_Population[1]->n++;
      /* The event occurs in the focal patch x: Recipient decreases */                   
      Positivity_Control(11, Table, x, 0, Y[R_0], J[R_0]);
      Y[R_0]--; J[R_0]--; Patch->n[0]--;                                                 
                         // Patch->Local_Strain_Population[0]->n--;
      
      Conjugation_Success = 1;                  
      x_Patch[2] = 0;   /* RECIPIENT */  x_Patch[3] = 1;   /* DONOR     */                    
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
  x_Patch[4] = Conjugation_Success;
  
  assert(Y[R_0] == (double)J[R_0] && J[R_0] == Patch->n[0]);
  assert(Y[R_1] == (double)J[R_1] && J[R_1] == Patch->n[1]);
  
  /* The event has been executed */
  
  /* (*Event) allows us to identify both the affected species ID 
     (donor, if the event has been a conjugation) and the patch where 
     the event has occurred, x */
  /* If the event is a between-patch movement, then x is the affected patch, 
     and y is the patch where the individual bacterial cell has moved to from patch x. 
     Otherwise, always x == y. 
  */
  /* If conjugation has occurred, the ID of recipient species is 0, 
     and the ID of the transconjugant species 1 
  */ 
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
  int check = 0;
  int n_Event; 
  int y, R, n_Patch;
  Community ** Patch = Table->Patch_System;

  int * J = Table->Vector_Model_Int_Variables;
  double * Y = Table->Vector_Model_Variables;

  if (Patch[x]->Out_Migration_Vector[Sp] != NULL && Patch[x]->No_NEI > 0) {
    n_Patch = Discrete_Sampling(Patch[x]->Out_Migration_Vector[Sp], Patch[x]->No_NEI) - 1;
  } 
  else {
    printf("Error: Invalid Out_Migration_Vector or No_NEI for patch %d and species %d\n", x, Sp);
    exit(1);
  }

  y = Patch[x]->Patch_Connections[n_Patch];
  R = y * Table->LOCAL_STATE_VARIABLES +  Sp; /* R is the model variable index 
                                                 of the species in the patch y */
  if (Sp == 0) 
    n_Event = 0; /* Outmigration (A_0) */
  else if (Sp == 1) 
    n_Event = 5; /* Outmigration (A_1) */
  else {
    printf("Species ID is ill-defined\n");
    printf("Sp = %d\n", Sp);
    printf("but Sp can be 1 or 0 in this context\n");
    printf("The program will exit\n");
    exit(1);
  }    
    
  check = Carrying_Capacity_Control_Check(n_Event, Table, y, Sp, Y[R], J[R]);
                  
  if (check == 1) {
    Y[R]++; J[R]++; Patch[x]->NEI[n_Patch]->n[Sp]++; Patch[x]->NEI[n_Patch]->m_0--;
     // Patch[x]->NEI[n_Patch]->Local_Strain_Population[Sp]->n++;
  }
  else {
    y = x; /* If the carrying capacity is reached, the individual dies in the transit
              to the other patch */
  }
   
  return (y);
}

