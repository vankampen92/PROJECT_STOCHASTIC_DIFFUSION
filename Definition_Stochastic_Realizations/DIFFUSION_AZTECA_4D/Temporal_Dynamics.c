/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2010 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

void Temporal_Dynamics(Community ** My_Community, Parameter_Table * Table, Stochastic_Rate * Rate)
{
  /* This function calculates the rates of all possible events from scratch, across and within cells 
     given certain system confituration as defined in My_Community 
  */
  int i,j,k,n, Sp;
  Community * P;
  int MODEL_STATE_VARIABLES;
  int No_of_CELLS;
  int No_of_EVENTS;
  int GRAND_No_of_EVENTS;
  double OutMigration;
  
  Parameter_Model * pa  = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
  
  P = My_Community[0];  /* P could be used as a pointer to the zero-th to be incremented 
			                     if necessary (not used like that in this implementation) 
		                    */
  No_of_CELLS             = pa->No_of_CELLS;
  MODEL_STATE_VARIABLES   = pa->MODEL_STATE_VARIABLES;
  Sp                      = pa->No_of_RESOURCES; 

  assert(Sp == 1);
  
  No_of_EVENTS            = pa->TOTAL_No_of_EVENTS; /* Total No of Events 
						                                           within each local population */  

  Immigration_Preassure_on_Focal_Patch_Initialization( My_Community, pa );

  Rate->max_Probability = 0.0;
  Rate->Total_Rate      = 0.0;

  /* K_W: Total Carrying Capacity (Workers), and Table->K_R is the max No of Worker per Nest */
  double K_W = (double)Table->K_R * (double)Table->Lambda_C_1;  /* -HK  -H7 */
  /* K_Q: Total Max No of Nests (per local patch) */ 
  double K_Q = (double)Table->Lambda_C_1;                       /* -H7      */

  GRAND_No_of_EVENTS = 0;
  for(i=0; i<No_of_CELLS; i++){

    P = My_Community[i];
    
    P->ratePatch = 0; 
    n = 0;
       
    /* 0: Worker Out-Migration (P --> P-1) and some other patch gains one */ 
    OutMigration = P->Total_Per_Capita_Out_Migration_Rate[W];
    assert( 4 * Table->Mu ==  OutMigration ); 
        
    P->rate[n] = OutMigration;                            P->rToI[n]  = OutMigration * (double)P->n[W]; 
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif 
    n++;

    /* 1: Worker External Immigration event */
    P->rate[n] = Table->Lambda_R_0;                        P->rToI[n]  = Table->Lambda_R_0; 
    P->ratePatch += P->rToI[n]; 
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;

    /* 2: Death of Workers  */
    P->rate[n] = Table->Delta_R_0;                         P->rToI[n]  = Table->Delta_R_0 * (double)P->n[W];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 3: Worker Production by Queens */
    P->rate[n] = Table->Beta_R * (K_W-(double)P->n[W])/K_W; P->rToI[n] = P->rate[n] * (double)P->n[Q];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 4: Nest Establishment  */
    P->rate[n] = Table->Eta_R *(K_Q-(double)P->n[Q])/K_Q;   P->rToI[n] = P->rate[n] * (double)P->n[W];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 5: Queen Death  */
    P->rate[n] = Table->Delta_R_1;                         P->rToI[n]  = Table->Delta_R_1 * (double)P->n[Q];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 6: Fly Out-Migration (F ---> F-1) and some other patch gains one */ 
    OutMigration = P->Total_Per_Capita_Out_Migration_Rate[F];
    assert( 4 * Table->Mu_C ==  OutMigration ); 
        
    P->rate[n] = OutMigration;                             P->rToI[n]  = OutMigration * (double)P->n[F];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 7: Fly External Immigration (F ---> F + 1) */  
    P->rate[n] = Table->Lambda_C_0;                        P->rToI[n]  = Table->Lambda_C_0; 
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 8: Fly Death  */
    P->rate[n] = Table->Delta_C_0;                         P->rToI[n]  = Table->Delta_C_0 * (double)P->n[F];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 9: Attack: Fly Consumption of a Worker Item and 'Dimmer' Formation */
    P->rate[n]= Table->Alpha_C_0*(double)P->n[W]/K_W;      P->rToI[n]= P->rate[n]*(double)P->n[F];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 10: Larval Development (within Worker heads) into adult flies  */ 
    P->rate[n] = Table->Nu_C_0;                            P->rToI[n] = P->rate[n]*(double)P->n[WF];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
    
    /* 11: Parasitized workers death */
    P->rate[n] = Table->Delta_C_1;                         P->rToI[n]  = Table->Delta_C_1 * (double)P->n[WF];
    P->ratePatch += P->rToI[n];
    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
    #endif
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
    #endif
    n++;
      
    #if defined BINARY_TREE_OPTIMIZATION
      Table->Leaves[i]->value = P->ratePatch;
    #endif

    assert( n == Table->TOTAL_No_of_EVENTS );
    
    Rate->Total_Rate += P->ratePatch;
    Rate->max_Probability = MAX( Rate->max_Probability, P->ratePatch );
  }

  #if defined BINARY_TREE_SUPER_OPTIMIZATION
  assert(GRAND_No_of_EVENTS == Table->TOTAL_GRAND_No_of_EVENTS);
  #endif 
  
  #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
  assert(GRAND_No_of_EVENTS == Table->TOTAL_GRAND_No_of_EVENTS);
  #endif 

  if(Rate->Total_Rate <= 0.0){
      printf("\n");
      printf(" R is the total temporal rate of system configuration change\n");
      printf(" R = %g\n", Rate->Total_Rate );
      printf(" As R is zero or negative, no change is possible\n");
      printf(" R shouldn't be negative. If it is, there are definitely some errors in the code\n");
      printf("\n");
      if( Rate->Total_Rate < 0.0 ) exit(0);
  }
}
