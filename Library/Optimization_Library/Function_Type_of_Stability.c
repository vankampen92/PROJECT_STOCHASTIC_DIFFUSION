#include <MODEL.h>

/* Under Constrution */

int Function_to_Type_of_Stability( Parameter_Table * Table )
{
  /* 	
   *  Input arguments:
   *
   *  . Table, a pointer to the main data structure controling all parameters
   *  of the excution.
  
   *  Output arguments:
   *
   *  . value, the value that takes this function for parameters values 
   *  as defined in input Table. 
   
   In this example, this is a three-valued function. Return values are

   . 0, for only resources or no resources at all. 

   . 1, for stable coexistce of resources and consumers

   . 2, for stable coexistce of resources and consumers and damped oscillations 
*/
  
  int i, j, k;  
  int Type_of_Stability; 

  Type_of_Stability = Coexistence_Condition ( Table ); 
  
  assert( Type_of_Stability )

  if (Type_of_Stability == 0 ) return(Type_of_Stability);
  else {
    
  }
  
  return(Type_of_Stability);
}


