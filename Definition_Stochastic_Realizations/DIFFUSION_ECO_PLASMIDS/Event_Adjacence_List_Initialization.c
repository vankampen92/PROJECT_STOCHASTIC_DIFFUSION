#include <MODEL.h>

void Event_Adjacence_List_Initialization(Community ** PATCH,
					 Parameter_Model * P)
{
  int i, j, m, no;
  
  no  = P->No_of_CELLS;
  m   = P->No_of_EVENTS; 

  int ** AD; 

  no  = P->No_of_CELLS;
 
  for(i=0; i<no; i++){
    
    AD = PATCH[i]->Event_Adjacence_List;
    for(j = 0; j < m; j++) {
      AD[0][j] = j;                               
      AD[1][j] = j;                               
      AD[2][j] = j;                               
      AD[3][j] = j;                               
      AD[4][j] = j; 
    }

    AD[0][m] = 7;   
    AD[1][m] = 7;
    AD[2][m] = 7;
    AD[3][m] = 7; 
    AD[4][m] = 7;                    
    AD[5][m] = 0;
    AD[6][m] = 0; 
  }
}



