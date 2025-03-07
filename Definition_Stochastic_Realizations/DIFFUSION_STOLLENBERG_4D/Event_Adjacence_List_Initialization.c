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

    AD[0][0] = 0; AD[0][1] = 2; AD[0][2] = 4;                                            AD[0][m] = 3;
    AD[1][0] = 0; AD[1][1] = 2; AD[1][2] = 4;                                            AD[1][m] = 3;
    AD[2][0] = 0; AD[2][1] = 2; AD[2][2] = 4;                                            AD[2][m] = 3;

    AD[3][0] = 0; AD[3][1] = 2; AD[3][2] = 4;                                            AD[3][m] = 3;
    AD[4][0] = 0; AD[4][1] = 2; AD[4][2] = 3; AD[4][3] = 4; AD[4][4] = 5; AD[4][5] = 9;  AD[4][m] = 6;
    AD[5][0] = 3; AD[5][1] = 4; AD[5][2] = 5; AD[5][3] = 9;                              AD[5][m] = 4;
    
    AD[6][0] = 6; AD[6][1] = 8; AD[6][2] = 9;                                            AD[6][m] = 3;
    AD[7][0] = 6; AD[7][1] = 8; AD[7][2] = 9;                                            AD[7][m] = 3;
    AD[8][0] = 6; AD[8][1] = 8; AD[8][2] = 9;                                            AD[8][m] = 3;

    AD[9][0] = 3; AD[9][1] = 4; AD[9][2] = 5; AD[9][3] = 6; AD[9][4] = 8; AD[9][5] = 9;
    AD[9][6] =10; AD[9][7] =11; AD[9][8] =12;                                            AD[9][m] = 9;

    AD[10][0]= 6; AD[10][1]= 8; AD[10][2]= 9; AD[10][3]=10; AD[10][4]=11; AD[10][5]=12;  AD[10][m]= 6;
    AD[11][0]= 6; AD[11][1]= 8; AD[11][2]= 9;                                            AD[11][m]= 3;
    AD[12][0]=10; AD[12][1]=11; AD[12][2]=12;                                            AD[12][m]= 3;

  }
}



