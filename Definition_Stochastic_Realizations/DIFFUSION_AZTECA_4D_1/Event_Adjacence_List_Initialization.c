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

    AD[0][0] = 0;  AD[0][1] = 2;  AD[0][2] = 7;  AD[0][3] = 11;                              AD[0][m] = 4;
    AD[1][0] = 0;  AD[1][1] = 2;  AD[1][2] = 7;  AD[1][3] = 11;                              AD[1][m] = 4;
    AD[2][0] = 0;  AD[2][1] = 2;  AD[2][2] = 7;  AD[2][3] = 11;                              AD[2][m] = 4;

    AD[3][0] = 3;  AD[3][1] = 5;  AD[3][2] = 6;  AD[3][3] = 7;                               AD[3][m] = 4;
    AD[4][0] = 3;  AD[4][1] = 5;  AD[4][2] = 6;  AD[4][3] = 7;                               AD[4][m] = 4;
    AD[5][0] = 3;  AD[5][1] = 5;  AD[5][2] = 6;  AD[5][3] = 7;                               AD[5][m] = 4;
    
    AD[6][0] = 0;  AD[6][1] = 2;  AD[6][2] = 7;  AD[6][3] = 11;                              AD[6][m] = 4;
    AD[7][0] = 3;  AD[7][1] = 5;  AD[7][2] = 6;  AD[7][3] = 7;                               AD[7][m] = 4;
  
    AD[8][0] = 8;  AD[8][1] = 10; AD[8][2] = 11;                                             AD[8][m] = 3;
    AD[9][0] = 8;  AD[9][1] = 10; AD[9][2] = 11;                                             AD[9][m] = 3;    
    AD[10][0]= 8;  AD[10][1]= 10; AD[10][2]= 11;                                             AD[10][m]= 3;
    
    AD[11][0]= 0;  AD[11][1]= 2;  AD[11][2]= 7;  AD[11][3]= 11; AD[11][4]= 12; AD[11][5]= 13; AD[11][m]= 6;
    AD[12][0]= 8;  AD[12][1]= 10; AD[12][2]= 11; AD[12][3]= 12; AD[12][4]= 13;                AD[12][m]= 5;
    AD[13][0]= 12; AD[13][1]= 13;                                                             AD[13][m]= 2;
  }
}



