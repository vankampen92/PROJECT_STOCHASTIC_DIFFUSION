#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH, Parameter_Model * P)
{
  /* 
     This is the subset of the Delta Matrix entries that are independent 
     from system configuration. Therefore, they only need to be 
     initialized once and for ever
  */
  
  int i, k, no, RP, A;
  double Out_Migration;

  Parameter_Table * Pa = (Parameter_Table *)P->Table; 

  double ** M;
  
  RP  = 0;    /* Resource Propagules */
  A   = 2;    /* Consumers           */
  no  = P->No_of_CELLS;
    
  for(i=0; i<no; i++){

    for(k=0; k < P->No_of_RESOURCES; k++){
 
      Out_Migration = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[k];
      assert ( PATCH[i]->No_NEI * P->Mu == Out_Migration );
    
      M = PATCH[i]->Event_Delta_Tensor[k];

      M[0][0] = -Out_Migration; M[0][1] = +Pa->Lambda_R_0; M[0][2] = -Pa->Delta_AP[k]; 
      M[1][0] = +Out_Migration; M[1][1] = -Pa->Lambda_R_0; M[1][2] = +Pa->Delta_AP[k];
      M[2][0] = -Out_Migration; M[2][1] = +Pa->Lambda_R_0; M[2][2] = -Pa->Delta_AP[k];
      M[3][0] = -Out_Migration; M[3][1] = +Pa->Lambda_R_0; M[3][2] = -Pa->Delta_AP[k];
      M[4][0] = +Out_Migration; M[4][1] = -Pa->Lambda_R_0; M[4][2] = +Pa->Delta_AP[k];   
    }
  }
}




