#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P)
{
  /* 
     This is the subset of the Delta Matrix entries that are independent 
     from system configuration. Therefore, they only need to be 
     initialized once and for ever
  */
  
  int i, no, RP, A;
  double Out_Migration_RP, Out_Migration_C;

  double ** M;
  
  RP  = 0;    /* Resource Propagules */
  A   = 2;    /* Consumers           */
  no  = P->No_of_CELLS;
    
  for(i=0; i<no; i++){
    Out_Migration_RP = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[RP];
    assert ( PATCH[i]->No_NEI * P->Mu == Out_Migration_RP );
    
    Out_Migration_C = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[A];
    assert( Out_Migration_C == PATCH[i]->No_NEI * P->Mu_C);
    
    M = PATCH[i]->Event_Delta_Matrix;

    M[0][0] = -Out_Migration_RP; M[0][2] = -P->Delta_R_1; 
    M[1][0] = +Out_Migration_RP; M[1][2] = +P->Delta_R_1;
    M[2][0] = -Out_Migration_RP; M[2][2] = -P->Delta_R_1;

    M[3][0] = +Out_Migration_RP; M[3][2] = +P->Delta_R_1;
    M[4][0] = -Out_Migration_RP; M[4][2] = -P->Delta_R_1;  M[4][3] = +P->Beta_R; M[4][5] = +P->Delta_R_0; 
    M[5][3] = -P->Beta_R;        M[5][5] = -P->Delta_R_0;

    M[6][6] = -Out_Migration_C;  M[6][8] = -P->Delta_C_0; 
    M[7][6] = +Out_Migration_C;  M[7][8] = +P->Delta_C_0; 
    M[8][6] = -Out_Migration_C;  M[8][8] = -P->Delta_C_0; 

    M[9][3] = -P->Beta_R;        M[9][5] = -P->Delta_R_0;
    M[9][6] = -Out_Migration_C;  M[9][8] = -P->Delta_C_0;
    M[9][10]= +P->Nu_C_0;        M[9][11]= +P->Beta_C;     M[9][12]= +P->Delta_C_0;
   
    M[10][6]= +Out_Migration_C;  M[10][8] = +P->Delta_C_0; 
    M[10][10]= -P->Nu_C_0;       M[10][11]= -P->Beta_C;    M[10][12]= -P->Delta_C_0;

    M[11][6]= +Out_Migration_C;  M[11][8] = +P->Delta_C_0;

    M[12][10]= -P->Nu_C_0;       M[12][11]= -P->Beta_C;    M[12][12]= -P->Delta_C_0;
  }
}




