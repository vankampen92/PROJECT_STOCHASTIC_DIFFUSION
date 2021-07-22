#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P)
{
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
    M[6][7] = +P->Nu_C_0;

    M[7][3] =+2.0*Out_Migration_C;                         M[7][5] =+2.0*P->Delta_C_0;
    M[7][7] =-P->Nu_C_0;

    M[8][0] = +Out_Migration_R; M[8][1] = -P->Lambda_R_0;  M[8][2] = +P->Delta_R_0;

    M[9][3] = -Out_Migration_C;                            M[9][5] = -P->Delta_C_0;
    M[9][7] = -P->Nu_C_0;                                  M[9][10]= +P->Eta_C_0;

    M[10][3] = +Out_Migration_C;                            M[10][5] = +P->Delta_C_0;
    M[10][7] = +P->Nu_C_0;                                 M[10][10]= -P->Eta_C_0;

  }
}




