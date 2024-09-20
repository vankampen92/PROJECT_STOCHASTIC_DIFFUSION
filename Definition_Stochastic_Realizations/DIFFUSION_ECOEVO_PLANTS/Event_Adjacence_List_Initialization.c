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

    AD[0][0] = 0; AD[0][1] = 2; AD[0][2] = 4; AD[0][3] = 5;                              AD[0][m] = 4;
    AD[1][0] = 0; AD[1][1] = 2; AD[1][2] = 4; AD[1][3] = 5;                              AD[1][m] = 4;
    AD[2][0] = 0; AD[2][1] = 2; AD[2][2] = 4; AD[2][3] = 5;                              AD[2][m] = 4;

    AD[3][0] = 0; AD[3][1] = 2; AD[3][2] = 4; AD[3][3] = 5;                              AD[3][m] = 4;

    AD[4][0] = 0; AD[4][1] = 2; AD[4][2] = 3; AD[4][3] = 4; AD[4][4] = 5; AD[4][5] = 7;  AD[4][m] = 6;
    AD[5][0] = 0; AD[5][1] = 2; AD[5][2] = 4; AD[5][3] = 5;                              AD[5][m] = 4;
    
    AD[6][0] = 3; AD[6][1] = 4; AD[6][2] = 5; AD[6][3] = 7;                              AD[6][m] = 4;

    AD[7][0] = 3; AD[7][1] = 4; AD[7][2] = 5; AD[7][3] = 7;                              AD[7][m] = 4;   
  }
}



