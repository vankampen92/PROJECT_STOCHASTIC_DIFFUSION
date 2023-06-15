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
  double OutMigration, n_R, n_RA, A_0, A_F, A_H, K_R;
  
  Parameter_Model * pa  = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
  
  P = My_Community[0]; /* P could be used as a pointer to the zero-th to be incremented 
			                    if necessary (not used like that in this implementation) 
		                   */
  No_of_CELLS             = pa->No_of_CELLS;
  MODEL_STATE_VARIABLES   = pa->MODEL_STATE_VARIABLES;
  Sp                      = pa->No_of_RESOURCES; 
  No_of_EVENTS            = pa->TOTAL_No_of_EVENTS; /* Total No of Events 
						                                           within each local population 
                                                    */  
  Immigration_Preassure_on_Focal_Patch_Initialization( My_Community, pa );

  Rate->max_Probability = 0.0;
  Rate->Total_Rate      = 0.0;

  K_R = (double)Table->K_R; 
  n_R = (double)Table->TOTAL_No_of_RESOURCES;
  A_0 = (double)Table->TOTAL_No_of_CONSUMERS;

  assert( No_of_CELLS == 1 );  /* The total number of consumers is globally constant, but
				                          the diffusion process changes this number locally. This 
				                          function calculates n_RA of the cell by considering 
				                          the total number of consumers is constant and equal 
				                          to A_0. This is only true if we have only one cell!!! 
			                          */
  for(i=0; i<No_of_CELLS; i++){
     
    P = My_Community[i];
    
    A_H = 0.0; 
    for(j=0; j<Sp; j++) A_H += (double)P->n[j]; 
    A_F = A_0 - A_H;    /* Number of Searching Consumers in a Current Patch */

    assert( A_F == (double)Table->TOTAL_No_of_FREE_CONSUMERS );

    P->ratePatch = 0; 
    n = 0;
    for(j=0; j<Sp; j++){
        
      /*   2*j: Attack: Consumer attack of a resource item and dimmer formation */
      P->rate[n]= Table->Theta_Consumers[j];  P->rToI[n]= P->rate[n]*A_F;
      P->ratePatch += P->rToI[n];
      n++;

      /* 2*j+1: Handling: Dimmer degrates and the individual consumer relaxes back 
                into searching mode) */
      P->rate[n] = Table->Nu_Consumers[j];    P->rToI[n]= P->rate[n]*(double)P->n[j]; 
      P->ratePatch += P->rToI[n];
      n++;
    }
      /*     n: Consumer Out-Migration (A --> A-1) and some other patch gains one */ 
      OutMigration = P->Total_Per_Capita_Out_Migration_Rate[0];
      assert( 4 * Table->Mu_C ==  OutMigration );
      assert( OutMigration == 0.0 );   /* No movement between patches in a one-patch system!!! */
        
      P->rate[n] = OutMigration;       P->rToI[n]  = OutMigration * A_F; 
      P->ratePatch += P->rToI[n];
      n++;

      /* n+1: Consumer External Immigration event  */
      assert(Table->Lambda_C_0 == 0.0);  /* No change of consumer population in 
                                            chemostatic conditions */  
      P->rate[n] = Table->Lambda_C_0;  P->rToI[n]  = Table->Lambda_C_0; 
      P->ratePatch += P->rToI[n];
      n++;
    
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
