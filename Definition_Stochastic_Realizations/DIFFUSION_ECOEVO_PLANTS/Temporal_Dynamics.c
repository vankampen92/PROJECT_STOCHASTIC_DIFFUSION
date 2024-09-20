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
  double K_R, m_0, y_S; 

  Parameter_Model * pa  = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites free from all species                   */
  
  P = My_Community[0];  /* P could be used as a pointer to the zero-th to be incremented 
			                     if necessary (not used like that in this implementation) 
		                    */
  No_of_CELLS             = pa->No_of_CELLS;
  MODEL_STATE_VARIABLES   = pa->MODEL_STATE_VARIABLES;
  Sp                      = pa->No_of_RESOURCES; 
  
  No_of_EVENTS            = pa->TOTAL_No_of_EVENTS; /* Total No of Events 
						                                           within each local population */  
  double * Y = Table->Vector_Model_Variables; 
  
  Immigration_Preassure_on_Focal_Patch_Initialization( My_Community, pa );

  Rate->max_Probability = 0.0;
  Rate->Total_Rate      = 0.0;
  
  GRAND_No_of_EVENTS = 0;
  for(i=0; i<No_of_CELLS; i++){

    P = My_Community[i];
    y_S = Local_Population_Resources(i, Y, Table);
    m_0 = K_R-y_S; 

    P->ratePatch = 0; 
    n = 0;
    for(k=0; k<Sp; k++) {

      RP = 2*k + Table->RP; 
      R  = 2*k + Table->R; 

      /* 0: Propagule Out-Migration (P --> P-1) and some other patch gains one */ 
      OutMigration = P->Total_Per_Capita_Out_Migration_Rate[RP];                    
      assert( 4 * Table->Mu ==  OutMigration ); 
          
      P->rate[n] = OutMigration;       P->rToI[n]  = OutMigration * (double)P->n[RP]; 
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif 
      n++;

      /* 1: Propagule External Immigration event */
      P->rate[n] = Table->Lambda_R_0;  P->rToI[n]  = Table->Lambda_R_0; 
      P->ratePatch += P->rToI[n]; 
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;

      /* 2: Propagule Death  */
      P->rate[n] = Table->Delta_R_1;   P->rToI[n]  = Table->Delta_RP[k] * (double)P->n[RP];
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
      
      /* 3: Propagule Production */
      P->rate[n] = Table->Beta_R;                                      P->rToI[n]  = Table->Beta_AP[k] * (double)P->n[R];
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
      
      /* 4: Propagule Establishment  (no mutation)*/
      P->rate[n] = Table->Eta_RP[k] * (1.0-2.0*Table->p_1) * m_0/K_R;  P->rToI[n]  = P->rate[n] * (double)P->n[RP];
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
      
      /* 5: Propagule Establishment (out mutation)  */
      P->rate[n] = Table->Eta_RP[k] * 2.0*Table->p_1 * m_0/K_R;        P->rToI[n]  = P->rate[n] * (double)P->n[RP];
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
      
      /* 6: Propagule Establishment (in mutation)   (always results from a mutation form another phenotype )*/  
      P->rate[n] = 0.0;                                                P->rToI[n]  = 0.0; 
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
      
      /* 7: Adult Death  A --> A - 1) */  
      P->rate[n] = Table->Delta_AP[k];                                 P->rToI[n]  = P->rate[n] * (double)P->n[R];
      P->ratePatch += P->rToI[n];
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
    }

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
