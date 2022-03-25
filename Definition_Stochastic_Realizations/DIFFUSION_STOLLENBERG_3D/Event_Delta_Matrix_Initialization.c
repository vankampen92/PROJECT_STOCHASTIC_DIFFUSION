#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P)
{
  /* 
     This is the subset of the Delta Matrix entries that are independent 
     from system configuration. Therefore, they only need to be 
     initialized 
  */
  
  int i, no, R, A;
  double Out_Migration_R, Out_Migration_C;

  double ** M; 
  R   = 0;    /* Resource */
  A   = 1;    /* Consumer */
  no  = P->No_of_CELLS;
    
  for(i=0; i<no; i++){
    Out_Migration_R = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[R];
    assert ( PATCH[i]->No_NEI * P->Mu == Out_Migration_R );
    
    Out_Migration_C = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[A];
    assert( Out_Migration_C == PATCH[i]->No_NEI * P->Mu_C);
    
    M = PATCH[i]->Event_Delta_Matrix;

    M[0][0] = -Out_Migration_R; M[0][1] = +P->Lambda_R_0;  M[0][2] = -P->Delta_R_0; 
    M[1][0] = +Out_Migration_R; M[1][1] = -P->Lambda_R_0;  M[1][2] = +P->Delta_R_0;
    M[2][0] = -Out_Migration_R; M[2][1] = +P->Lambda_R_0;  M[2][2] = -P->Delta_R_0;

    M[3][3] = -Out_Migration_C;                            M[3][5] = -P->Delta_C_0; 
    M[4][3] = +Out_Migration_C;                            M[4][5] = +P->Delta_C_0;
    M[5][3] = -Out_Migration_C;                            M[5][5] = -P->Delta_C_0;

    M[6][0] = -Out_Migration_R; M[6][1] = +P->Lambda_R_0;  M[6][2] = -P->Delta_R_0; 
    M[6][3] = -Out_Migration_C;                            M[6][5] = -P->Delta_C_0;
    M[6][8] = +P->Beta_C;       M[6][9] = +P->Delta_C_0;   M[6][10] = +P->Nu_C_0;

    M[7][0] = +Out_Migration_R; M[7][1] = -P->Lambda_R_0;  M[7][2] = +P->Delta_R_0;

    M[8][3] = +Out_Migration_C;                            M[8][5] = +P->Delta_C_0;

    M[9][8] = -P->Beta_C;       M[9][9] = -P->Delta_C_0;   M[6][10]= -P->Nu_C_0;

    M[10][3] = +Out_Migration_C;                           M[10][5] = +P->Delta_C_0;
    M[10][8] = -P->Beta_C;      M[10][9] = -P->Delta_C_0;  M[10][10]= -P->Nu_C_0;
  }
}




