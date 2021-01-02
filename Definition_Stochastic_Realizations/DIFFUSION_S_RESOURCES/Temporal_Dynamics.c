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
  double OutMigration;
  
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
						       within each local population */  

  Immigration_Preassure_on_Focal_Patch_Initialization( My_Community, pa );

  Rate->max_Probability = 0.0;
  Rate->Total_Rate      = 0.0;
  
  for(i=0; i<No_of_CELLS; i++){

    P = My_Community[i];
    
    P->ratePatch = 0; 
    n = 0; 
    for(j=0; j<Sp; j++) {
    
      /* 0: Out-Migration (S --> S-1) and some other patch gains one */ 
      OutMigration = P->Total_Per_Capita_Out_Migration_Rate[j];

      assert( 4 * Table->Mu ==  OutMigration ); 
        
      /* Probability rate for each of the events */     
      P->rate[n] = OutMigration;     P->rToI[n]  = OutMigration * (double)P->n[j];
      
      P->ratePatch += P->rToI[n];

      n++;
      
      /* 1 : External Immigration event of species j */
      
      P->rate[n] = Table->Lambda_R[j];  P->rToI[n]  = Table->Lambda_R[j]; 
      
      P->ratePatch += P->rToI[n];

      n++;

      /* 2 : Exponential Decaly of Resource j */
      
      P->rate[n] = Table->Delta_R[j];  P->rToI[n]  = Table->Delta_R[j] * (double)P->n[j];
      
      P->ratePatch += P->rToI[n];

      n++;
    }    

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
