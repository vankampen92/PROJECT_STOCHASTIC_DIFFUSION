#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P)
{
  /* 
     This is the subset of the Delta Matrix entries that are independent 
     from system configuration. Therefore, they only need to be 
     initialized once and for ever
  */
  
  int i, no, W, Q, F;
  double Out_Migration_W, Out_Migration_Q, Out_Migration_F;

  double ** M;
  
  W  = 0;    /* Workers */
  Q  = 1;    /* Queens  */
  F  = 2;    /* Flies   */
  no  = P->No_of_CELLS;
    
  for(i=0; i<no; i++){
    Out_Migration_W = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[W];
    assert (Out_Migration_W == PATCH[i]->No_NEI * P->Mu);

    Out_Migration_Q = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[Q];
    assert( Out_Migration_Q == PATCH[i]->No_NEI * P->Lambda_C_1);
    
    Out_Migration_F = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[F];
    assert( Out_Migration_F == PATCH[i]->No_NEI * P->Mu_C);
    
    M = PATCH[i]->Event_Delta_Matrix;

    M[0][0] = -Out_Migration_W; M[0][2]  = -P->Delta_R_0; 
    M[1][0] = +Out_Migration_W; M[1][2]  = +P->Delta_R_0;
    M[2][0] = -Out_Migration_W; M[2][2]  = -P->Delta_R_0;

    M[3][3] = -Out_Migration_Q; M[3][5]  = -P->Delta_R_1; M[3][6] = -P->Beta_R;
    M[4][3] = +Out_Migration_Q; M[4][5]  = +P->Delta_R_1; M[4][6] = +P->Beta_R;
    M[5][3] = -Out_Migration_Q; M[5][5]  = -P->Delta_R_1; M[5][6] = -P->Beta_R; 

    M[6][0] = +Out_Migration_W; M[6][2]  = +P->Delta_R_0;  
    M[7][3] = +Out_Migration_Q; M[7][5]  = +P->Delta_R_1; M[7][6] = +P->Beta_R; 

    M[8][8]  = -Out_Migration_F;  M[8][10]  = -P->Delta_C_0; 
    M[9][8]  = +Out_Migration_F;  M[9][10]  = +P->Delta_C_0; 
    M[10][8] = -Out_Migration_F;  M[10][10] = -P->Delta_C_0; 

    M[11][0] = -Out_Migration_W;  M[11][2]  = -P->Delta_R_0; M[11][12] = +P->Nu_C_0; M[11][13]= +P->Delta_C_1;
    M[12][8] = +Out_Migration_F;  M[12][10] = +P->Delta_C_0; M[12][12] = -P->Nu_C_0; M[12][13]= -P->Delta_C_1;    
    M[13][12]= -P->Nu_C_0;        M[13][13] = -P->Delta_C_1;
  }
}




