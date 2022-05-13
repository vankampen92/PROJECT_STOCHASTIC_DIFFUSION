#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH,
				       Parameter_Model * P)
{
  int i, no, R, A;
  double Out_Migration_C;
  double K_R, n_R; 

  double ** M; 

  A   = 0;    /* Consumer */ /* Resources are fixed to a density n_R/K_R */
  no  = P->No_of_CELLS;


  K_R = (double)P->K_R; 
  n_R = (double)P->TOTAL_No_of_RESOURCES; 
    
  for(i=0; i<no; i++){
   
    Out_Migration_C = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[A];
    assert( Out_Migration_C == PATCH[i]->No_NEI * P->Mu_C);
    
    M = PATCH[i]->Event_Delta_Matrix;

    M[0][0] = -Out_Migration_C; M[0][2] = -P->Alpha_C_0 * n_R/K_R; 


    M[1][0] = +Out_Migration_C; M[1][2] = +P->Alpha_C_0 * n_R/K_R;
    

    M[2][0] = -Out_Migration_C; M[2][2] = -P->Alpha_C_0 * n_R/K_R;  M[2][3] = +P->Nu_C_0;
    

    M[3][0] = +Out_Migration_C; M[3][2] = +P->Alpha_C_0 * n_R/K_R;  M[3][3] = -P->Nu_C_0;


    M[4][0] = -Out_Migration_C; M[4][2] = -P->Alpha_C_0 * n_R/K_R;  M[4][3] = -P->Nu_C_0;  M[4][5] = P->Eta_C_0;
    

    M[5][0] = +Out_Migration_C; M[5][2] = +P->Alpha_C_0 * n_R/K_R;  M[5][3] = +P->Nu_C_0;  M[5][5] = -P->Eta_C_0;  
    
  }
}




