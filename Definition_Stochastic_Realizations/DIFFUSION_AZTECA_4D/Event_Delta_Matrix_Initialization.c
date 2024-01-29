#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P)
{
  /* 
     This is the subset of the Delta Matrix entries that are independent 
     from system configuration. Therefore, they only need to be 
     initialized once and for ever
  */
  
  int i, no, W, F;
  double Out_Migration_W, Out_Migration_F;

  double ** M;
  
  W  = 0;    /* Workers */
  F  = 2;    /* Flies   */
  no  = P->No_of_CELLS;
    
  for(i=0; i<no; i++){
    Out_Migration_W = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[W];
    assert ( PATCH[i]->No_NEI * P->Mu == Out_Migration_W );
    
    Out_Migration_F = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[F];
    assert( Out_Migration_F == PATCH[i]->No_NEI * P->Mu_C);
    
    M = PATCH[i]->Event_Delta_Matrix;

    M[0][0] = -Out_Migration_W; M[0][2]  = -P->Delta_R_0; 
    M[1][0] = +Out_Migration_W; M[1][2]  = +P->Delta_R_0;
    M[2][0] = -Out_Migration_W; M[2][2]  = -P->Delta_R_0;

    M[3][0] = +Out_Migration_W; M[3][2]  = +P->Delta_R_0;
    M[4][0] = -Out_Migration_W; M[4][2]  = -P->Delta_R_0;  M[4][5] = +P->Delta_R_1; 
    M[5][5] = -P->Delta_R_1;

    M[6][6] = -Out_Migration_F;  M[6][8]  = -P->Delta_C_0; 
    M[7][6] = +Out_Migration_F;  M[7][8]  = +P->Delta_C_0; 
    M[8][6] = -Out_Migration_F;  M[8][8]  = -P->Delta_C_0; 

    M[9][0] = -Out_Migration_W;  M[9][2]  = -P->Delta_R_0; M[9][10] = +P->Nu_C_0;  M[9][11]= +P->Delta_C_1;
    M[10][6]= +Out_Migration_F;  M[10][8] = +P->Delta_C_0; M[10][10]= -P->Nu_C_0; M[10][11]= -P->Delta_C_1;    
    M[11][10]= -P->Nu_C_0;       M[11][11]= -P->Delta_C_1;
  }
}




