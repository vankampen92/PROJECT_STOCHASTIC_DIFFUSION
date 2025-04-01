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
*/
  
int Type_of_Stability; 

    Type_of_Stability = (int)Calculate_Stability_Stationary_Point( Table );
    
return(Type_of_Stability);
}

double Function_to_Type_of_Stability_Double( Parameter_Table * Table )
{
  /* 	
   *  Input arguments:
   *
   *  . Table, a pointer to the main data structure controling all parameters
   *  of the excution.
  
   *  Output arguments:
   *
   *  . Type_of_Stability, the value that takes this function for parameters values 
   *  as defined in input Table. 
*/  
  double Type_of_Stability; 
    
  Type_of_Stability = Calculate_Stability_Stationary_Point( Table );    
  
return(Type_of_Stability);
}

