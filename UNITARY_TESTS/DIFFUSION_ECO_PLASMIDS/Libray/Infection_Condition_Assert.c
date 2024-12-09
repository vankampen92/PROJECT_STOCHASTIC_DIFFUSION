#include <MODEL.h>

int Infection_Condition_Assert(int i_Sp, int * Profile, Parameter_Table * Table)
{
  int k; 
  int n_bool; 
  int n_count; 
  int S; 

  S = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++) 
      if( Profile[k] == 1) S++;

    /* Counting infections */
  n_count = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++) {
      if( Profile[k] == 1) {  
            /* Infection Constraint */ 
            if(Table->IBP[i_Sp][k] == 1.0) n_count++; 
      }
    }
    
  if(n_count == S) n_bool = 1;

  return(n_bool); 
}
