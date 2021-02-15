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

    AD[0][0] = 0; AD[0][1] = 1; AD[0][2] = 2; AD[0][3] = 6; AD[0][4] = 8;                AD[0][m] = 5;
    AD[1][0] = 0; AD[1][1] = 1; AD[1][2] = 2; AD[1][3] = 6; AD[1][4] = 8;                AD[1][m] = 5;
    AD[2][0] = 0; AD[2][1] = 1; AD[2][2] = 2; AD[2][3] = 6; AD[2][4] = 8;                AD[2][m] = 5;

    AD[3][0] = 3; AD[3][1] = 5; AD[3][2] = 6; AD[3][3] = 9;                              AD[3][m] = 4;
    AD[4][0] = 3; AD[4][1] = 5; AD[4][2] = 6; AD[4][3] = 9;                              AD[4][m] = 4;
    AD[5][0] = 3; AD[5][1] = 5; AD[5][2] = 6; AD[5][3] = 9;                              AD[5][m] = 4;

    AD[6][0] = 0; AD[6][1] = 1; AD[6][2] = 2; AD[6][3] = 3; AD[6][4] = 5; AD[6][5] = 6;
    AD[6][6] = 7; AD[6][7] = 8; AD[6][8] = 9;                                            AD[6][m] = 9;

    AD[7][0] = 3; AD[7][1] = 5; AD[7][2] = 6; AD[7][3] = 7; AD[7][4] = 9;                AD[7][m] = 5;
    AD[8][0] = 0; AD[8][1] = 1; AD[8][2] = 2; AD[8][3] = 6; AD[8][4] = 8;                AD[8][m] = 5;
    AD[9][0] = 3; AD[9][1] = 5; AD[9][2] = 6; AD[9][3] = 7; AD[9][4] = 9; AD[9][5] = 10; AD[9][m] = 6;

    AD[10][0]= 3; AD[10][1]= 5; AD[10][2]= 6; AD[10][3]= 7; AD[10][4]= 9; AD[10][5]= 10; AD[10][m]= 6;
  }
}



