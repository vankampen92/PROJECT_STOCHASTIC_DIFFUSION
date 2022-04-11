#include <MODEL.h>

void Event_Adjacence_List_Initialization(Community ** PATCH,
					 Parameter_Model * P)
{
  int i, m, no;
  
  no  = P->No_of_CELLS;
  m   = P->TOTAL_No_of_EVENTS; 

  int ** AD; 

  no  = P->No_of_CELLS;
 
  for(i=0; i<no; i++){
    
    AD = PATCH[i]->Event_Adjacence_List;

    AD[0][0] = 0; AD[0][1] = 2;                 AD[0][m] = 2;

    AD[1][0] = 0; AD[1][1] = 2;                 AD[1][m] = 2;

    AD[2][0] = 0; AD[2][1] = 2; AD[2][2] = 3;   AD[2][m] = 3;

    AD[3][0] = 0; AD[3][1] = 2; AD[3][2] = 3;   AD[3][m] = 3;

  }
}



