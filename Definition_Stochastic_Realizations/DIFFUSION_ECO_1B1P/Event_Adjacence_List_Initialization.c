#include <MODEL.h>

void Event_Adjacence_List_Initialization(Community ** PATCH, Parameter_Model * P)
{
  int i, j, m, no;
  int ** AD;
  
  m   = P->TOTAL_No_of_EVENTS; 
    
  no  = P->No_of_CELLS;

  for(i=0; i<no; i++){
    
    AD = PATCH[i]->Event_Adjacence_List;

    for(j = 0; j < m; j++) 
      AD[j][m] = 10; /* All events change system configuration in a way that there are a total of 10 other event rates affected */

    AD[0][0] = 0; AD[0][1] = 1; AD[0][2] = 2; AD[0][3] = 3; AD[0][4] = 4; AD[0][5] = 6; AD[0][6] = 8; AD[0][7] = 9; AD[0][8] = 10; AD[0][9] = 11;
    AD[2][0] = 0; AD[2][1] = 1; AD[2][2] = 2; AD[2][3] = 3; AD[2][4] = 4; AD[2][5] = 6; AD[2][6] = 8; AD[2][7] = 9; AD[2][8] = 10; AD[2][9] = 11;
    AD[3][0] = 0; AD[3][1] = 1; AD[3][2] = 2; AD[3][3] = 3; AD[3][4] = 4; AD[3][5] = 6; AD[3][6] = 8; AD[3][7] = 9; AD[3][8] = 10; AD[3][9] = 11;

    AD[5][0] = 1; AD[5][1] = 3; AD[5][2] = 4; AD[5][3] = 5; AD[5][4] = 6; AD[5][5] = 7; AD[5][6] = 8; AD[5][7] = 9; AD[5][8] = 10; AD[5][9] = 11;
    AD[7][0] = 1; AD[7][1] = 3; AD[7][2] = 4; AD[7][3] = 5; AD[7][4] = 6; AD[7][5] = 7; AD[7][6] = 8; AD[7][7] = 9; AD[7][8] = 10; AD[7][9] = 11;
    AD[8][0] = 1; AD[8][1] = 3; AD[8][2] = 4; AD[8][3] = 5; AD[8][4] = 6; AD[8][5] = 7; AD[8][6] = 8; AD[8][7] = 9; AD[8][8] = 10; AD[8][9] = 11;
        
    AD[6][0] = 1; AD[6][1] = 3; AD[6][2] = 4; AD[6][3] = 5; AD[6][4] = 6; AD[6][5] = 7; AD[6][6] = 8; AD[6][7] = 9; AD[6][8] = 10; AD[6][9] = 11;
    AD[9][0] = 1; AD[9][1] = 3; AD[9][2] = 4; AD[9][3] = 5; AD[9][4] = 6; AD[9][5] = 7; AD[9][6] = 8; AD[9][7] = 9; AD[9][8] = 10; AD[9][9] = 11;

    AD[11][0] = 0; AD[11][1] = 2; AD[11][2] = 3; AD[11][3] = 4; AD[11][4] = 5; AD[11][5] = 7; AD[11][6] = 8; AD[11][7] = 9; AD[11][8] = 10; AD[11][9] = 11;

     AD[1][0] = 0;  AD[1][1] = 1;  AD[1][2] = 2;  AD[1][3] = 3;  AD[1][4] = 4;  AD[1][5] = 6;  AD[1][6] = 8;  AD[1][7] = 9;  AD[1][8] = 10;  AD[1][9] = 11;
     AD[4][0] = 0;  AD[4][1] = 1;  AD[4][2] = 2;  AD[4][3] = 3;  AD[4][4] = 4;  AD[4][5] = 6;  AD[4][6] = 8;  AD[4][7] = 9;  AD[4][8] = 10;  AD[4][9] = 11;
    AD[10][0] = 0; AD[10][1] = 1; AD[10][2] = 2; AD[10][3] = 3; AD[10][4] = 4; AD[10][5] = 6; AD[10][6] = 8; AD[10][7] = 9; AD[10][8] = 10; AD[10][9] = 11;
  }
}



