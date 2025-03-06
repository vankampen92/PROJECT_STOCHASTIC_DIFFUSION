#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH, Parameter_Model * P)
{
  /* 
     This is the subset of the Delta Matrix entries that are independent 
     from system configuration. Therefore, they only need to be 
     initialized once and for ever
  */
  
  int i, no;
  double Eps;
  double d_0, d_1;
  double mu_0, mu_1;
  double la_0, la_1; 
  double ** M;
  
  d_0 = P->Delta_AP[0];
  d_1 = P->Delta_AP[1];
  
  la_0 = P->Lambda_R_0; /* Both strains receive external immigrants at the same rate */
  la_1 = P->Lambda_R_0;
  Eps = P->p_1;
  
  no  = P->No_of_CELLS;

  for(i=0; i<no; i++){
    mu_0 = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[0];
    assert( PATCH[i]->No_NEI * P->Mu_RP[0] == mu_0 );
    
    mu_1 = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[1];
    assert( PATCH[i]->No_NEI * P->Mu_RP[1] == mu_1 );

    M = PATCH[i]->Event_Delta_Matrix;

    M[0][0] = -mu_0; M[0][1] = +la_0; M[0][2] = -d_0; M[0][6] = +la_1;
    M[2][0] = -mu_0; M[2][1] = +la_0; M[2][2] = -d_0; M[2][6] = +la_1;
    M[3][0] = -mu_0; M[3][1] = +la_0; M[3][2] = -d_0; M[3][6] = +la_1;

    M[1][0] = +mu_0; M[1][1] = -la_0; M[1][2] = d_0; M[1][6] = -la_1;
    M[4][0] = +mu_0; M[4][1] = -la_0; M[4][2] = d_0; M[4][6] = -la_1;
    M[10][0]= +mu_0; M[10][1]= -la_0; M[10][2]= d_0; M[10][6]= -la_1;

    M[5][1] = +la_0; M[5][5] = -mu_1; M[5][6] = +la_1; M[5][7] = -d_1;
    M[7][1] = +la_0; M[7][5] = -mu_1; M[7][6] = +la_1; M[7][7] = -d_1;
    M[8][1] = +la_0; M[8][5] = -mu_1; M[8][6] = +la_1; M[8][7] = -d_1;

    M[6][1] = -la_0; M[6][5] = +mu_1; M[6][6] = -la_1; M[6][7] = +d_1;
    M[9][1] = -la_0; M[9][5] = +mu_1; M[9][6] = -la_1; M[9][7] = +d_1;

    M[11][0]= -mu_0; M[11][2] = -d_0; M[11][5]= +mu_1; M[11][7] = +d_1;   
  }
}




