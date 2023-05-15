/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2010 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

void Temporal_Dynamics(Community ** My_Community, Parameter_Table * Table, Stochastic_Rate * Rate)
{
  int i,j,k,n, Sp;
  Community * P;
  int MODEL_STATE_VARIABLES;
  int No_of_CELLS;
  int No_of_EVENTS;
  double OutMigration, n_R, n_RA, A_0, K_R;
  
  Parameter_Model * pa  = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
  
  P = My_Community[0]; /* P could be used as a pointer to the zero-th to be incremented 
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
    
    P->ratePatch = 0; 
    n = 0;

    n_RA = (A_0 - (double)P->n[A]);  
    
    /* 0: Consumer Out-Migration (A --> A-1) and some other patch gains one */ 
    OutMigration = P->Total_Per_Capita_Out_Migration_Rate[A];
    assert( 4 * Table->Mu_C ==  OutMigration ); 
        
    P->rate[n] = OutMigration;       P->rToI[n]  = OutMigration * (double)P->n[A];
    P->ratePatch += P->rToI[n];
    n++;

    /* 1: Consumer External Immigration event  */  
    P->rate[n] = Table->Lambda_C_0;  P->rToI[n]  = Table->Lambda_C_0; 
    P->ratePatch += P->rToI[n];
    n++;

    /* 2: Consumer Consumption of Resource and dimmer formation */
    P->rate[n]= Table->Alpha_C_0* n_R/K_R;     P->rToI[n]= P->rate[n]*(double)P->n[A];
    P->ratePatch += P->rToI[n];
    n++;

    /* 3: Dimmer degrations (production of two consumers */
    P->rate[n] = Table->Nu_C_0;                           P->rToI[n]= P->rate[n]*n_RA; 
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
