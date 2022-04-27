#include <MODEL.h>

void Common_Initial_Condition_Command_Line_Arguments_into_Table(Parameter_Table * Table)
{
    /*  This definition is contingent to TYPE of MODEL at work from the pre-defined family of models:
	DIFFUSION_BD_2D, DIFFUSION_HII_1D, ...

	This definitions can be taken by both deterministic functions (associated to the ODE system
	numerical integration) or the stochastic functions (related to either the generation of 
	stochastic replicates of the integration of the master equation). 
    */ 
    Table->TOTAL_No_of_RESOURCES  = (int)(Table->p_1 * (double)Table->K_R);
    Table->TOTAL_No_of_CONSUMERS  = Table->No_of_INDIVIDUALS;  /* -HN 20 as input argument */ 
    
    assert(Table->p_2 <= 1.0 && Table->p_2 >= 0.0);  // 
    assert(Table->p_1 <= 1.0 && Table->p_1 >= 0.0);  // Fractions!!!  
    
    Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0 = (int)(Table->p_2*(double)Table->TOTAL_No_of_CONSUMERS);
    Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = Table->TOTAL_No_of_CONSUMERS-Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
}
  
void Initial_Condition_Master_Equation( Parameter_Table * Table, double * y_INI )
{
  int i, i_0, No_of_CONFIGURATIONAL_STATES; 
  int n_Time_0, m_Time_0;
  int a_0;
  int n, m;

  No_of_CONFIGURATIONAL_STATES = Table->MEq->No_of_CONFIGURATIONAL_STATES;
  a_0   = Table->TOTAL_No_of_CONSUMERS;

  /* No of Free Consumers at time 0          */         /* n */ 
  n_Time_0 = Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  /* No of Handling Consumers at time 0      */         /* m */
  m_Time_0 = a_0 - Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;

  /* This initial condition involves no triplets at time t = 0.0 */
  
  assert(m_Time_0 + n_Time_0 <= a_0); 

  nm_to_i_Map(Table,
	      &i_0, n_Time_0, m_Time_0);  
 
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {
    
    if( i == i_0 ) y_INI[i] = 1.0;
    else           y_INI[i] = 0.0;       
    
  }
}

							 
  
  
  
