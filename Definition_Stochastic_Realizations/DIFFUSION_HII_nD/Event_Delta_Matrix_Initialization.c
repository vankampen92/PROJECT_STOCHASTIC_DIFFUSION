#include <MODEL.h>

void Event_Delta_Matrix_Initialization(Community ** PATCH, Parameter_Model * P)
{
  int i, j, k, m, no, Sp, R;
  int TOTAL_No_of_EVENTS; 
  double Out_Migration_C;
  double K_R, n_R; 

  double ** M; 

  no  = P->No_of_CELLS;
  m   = P->TOTAL_No_of_EVENTS; 
  Sp  = P->No_of_RESOURCES; 

  assert(no == 1);
  
  for(i=0; i<no; i++){
   
    M = PATCH[i]->Event_Delta_Matrix;

    TOTAL_No_of_EVENTS = 0; 
    for(j=0; j<Sp; j++) {

      Out_Migration_C = PATCH[i]->Total_Per_Capita_Out_Migration_Rate[j];
      assert( Out_Migration_C == PATCH[i]->No_NEI * P->Mu_C );
      assert( Out_Migration_C == 0.0 ); /* No Out-Migration (movement) of handling consumers 
                                           from j = 0, ..., Sp-1
                                        */
      /* Attack   */
      for(k=0; k<Sp; k++) { 
        M[j*P->No_of_EVENTS][k*P->No_of_EVENTS]   = -P->Theta_Consumers[k];
      }
      M[j*P->No_of_EVENTS][j*P->No_of_EVENTS + 1] = +P->Nu_Consumers[j];
      M[j*P->No_of_EVENTS][Sp*P->No_of_EVENTS]    = -Out_Migration_C;
      TOTAL_No_of_EVENTS++;

      /* Handling */
      for(k=0; k<Sp; k++) { 
        M[j*P->No_of_EVENTS+1][k*P->No_of_EVENTS]   = +P->Theta_Consumers[k];
      }
      M[j*P->No_of_EVENTS+1][j*P->No_of_EVENTS + 1] = -P->Nu_Consumers[j];
      M[j*P->No_of_EVENTS+1][Sp*P->No_of_EVENTS]    = +Out_Migration_C;  
      TOTAL_No_of_EVENTS++;  
    }   

    /* Out-Migration */
    for(k=0; k<Sp; k++) { 
        M[Sp*P->No_of_EVENTS][k*P->No_of_EVENTS]  = -P->Theta_Consumers[k];
    }
    M[Sp*P->No_of_EVENTS][Sp*P->No_of_EVENTS]     = -Out_Migration_C;
    TOTAL_No_of_EVENTS++;

    /* In-Migration */
    for(k=0; k<Sp; k++) { 
        M[Sp*P->No_of_EVENTS][k*P->No_of_EVENTS]  = +P->Theta_Consumers[k];
    }
    M[Sp*P->No_of_EVENTS][Sp*P->No_of_EVENTS]     = +Out_Migration_C;
    TOTAL_No_of_EVENTS++;

    assert(TOTAL_No_of_EVENTS == m);
  }
}




