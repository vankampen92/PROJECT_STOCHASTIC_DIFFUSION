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

  double K_R = (double)Table->K_R; 
  
  for(i=0; i<No_of_CELLS; i++){

    P = My_Community[i];
    
    P->ratePatch = 0; 
    n = 0;
       
    /* 0: Propagule Out-Migration (P --> P-1) and some other patch gains one */ 
    OutMigration = P->Total_Per_Capita_Out_Migration_Rate[RP];
    assert( 4 * Table->Mu ==  OutMigration ); 
        
    P->rate[n] = OutMigration;       P->rToI[n]  = OutMigration * (double)P->n[RP]; 
    P->ratePatch += P->rToI[n];
    n++;
      
    /* 1: Propagule External Immigration event */
    P->rate[n] = Table->Lambda_R_0;  P->rToI[n]  = Table->Lambda_R_0; 
    P->ratePatch += P->rToI[n];
    n++;
    
    /* 2: Propagule Death  */
    P->rate[n] = Table->Delta_R_1;   P->rToI[n]  = Table->Delta_R_1 * (double)P->n[RP];
    P->ratePatch += P->rToI[n];
    n++;

    /* 3: Propagule Production */
    P->rate[n] = Table->Beta_R;                            P->rToI[n]  = Table->Beta_R * (double)P->n[R];
    P->ratePatch += P->rToI[n];
    n++;

    /* 4: Propagule Establishment  */
    P->rate[n] = Table->Eta_R *(K_R-(double)P->n[R])/K_R;  P->rToI[n]  = P->rate[n] * (double)P->n[RP];
    P->ratePatch += P->rToI[n];
    n++;

    /* 5: Resource Death  */
    P->rate[n] = Table->Delta_R_0;                         P->rToI[n]  = Table->Delta_R_0 * (double)P->n[R];
    P->ratePatch += P->rToI[n];
    n++;
    
    /* 6: Consumer Out-Migration (A ---> A-1) and some other patch gains one */ 
    OutMigration = P->Total_Per_Capita_Out_Migration_Rate[A];
    assert( 4 * Table->Mu_C ==  OutMigration ); 
        
    P->rate[n] = OutMigration;                             P->rToI[n]  = OutMigration * (double)P->n[A];
    P->ratePatch += P->rToI[n];
    n++;
      
    /* 7: Consumer External Immigration (A ---> A + 1) */  
    P->rate[n] = Table->Lambda_C_0;                        P->rToI[n]  = Table->Lambda_C_0; 
    P->ratePatch += P->rToI[n];
    n++;

    /* 8: Searching Consumer Death  */
    P->rate[n] = Table->Delta_C_0;                         P->rToI[n]  = Table->Delta_C_0 * (double)P->n[A];
    P->ratePatch += P->rToI[n];
    n++;
    
    /* 9: Attack: Consumer Consumption of a Resource Item and Dimmer Formation */
    P->rate[n]= Table->Alpha_C_0*(double)P->n[R]/K_R;      P->rToI[n]= P->rate[n]*(double)P->n[A];
    P->ratePatch += P->rToI[n];
    n++;

    /* 10: Handling Consumers relax back into Free Consumers */ 
    P->rate[n] = Table->Nu_C_0;                            P->rToI[n] = P->rate[n]*(double)P->n[RA];
    P->ratePatch += P->rToI[n];
    n++;

    /* 11: Production of new Consumer Individuals */
    P->rate[n] = Table->Beta_C;                            P->rToI[n] = P->rate[n]*(double)P->n[RA]; 
    P->ratePatch += P->rToI[n];
    n++;

    /* 12: Handling Consumer Death */
    P->rate[n] = Table->Delta_C_0;                         P->rToI[n]  = Table->Delta_C_0 * (double)P->n[RA];
    P->ratePatch += P->rToI[n];
    n++;

    #if defined BINARY_TREE_OPTIMIZATION
    Table->Leaves[i]->value = P->ratePatch;
    #endif

    assert( n == Table->TOTAL_No_of_EVENTS );
    
    Rate->Total_Rate += P->ratePatch;
    Rate->max_Probability = MAX( Rate->max_Probability, P->ratePatch );
  }

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
