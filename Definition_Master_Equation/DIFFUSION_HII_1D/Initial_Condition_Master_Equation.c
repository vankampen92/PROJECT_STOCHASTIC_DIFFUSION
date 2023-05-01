#include <MODEL.h>

void Initial_Condition_Master_Equation( Parameter_Table * Table, double * y_INI )
{
  int i, No_of_CONFIGURATIONAL_STATES; 
  int m_Time_0;
  int n_0, a_0;
  int n;

  No_of_CONFIGURATIONAL_STATES = Table->MEq->No_of_CONFIGURATIONAL_STATES;
  a_0   = Table->TOTAL_No_of_CONSUMERS;

  /* No of Free Consumers at time 0          */
  m_Time_0 = Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0;
  /* No accumulated Feeding Events at time 0 */
  n_0      = 0;                                        
  
  assert(m_Time_0 <= a_0); 
 
  for (i=0; i<No_of_CONFIGURATIONAL_STATES; i++) {
    
    if( i == m_Time_0 ) y_INI[i] = 1.0;
    else                y_INI[i] = 0.0;       
    
  }
}

							 
  
  
  
