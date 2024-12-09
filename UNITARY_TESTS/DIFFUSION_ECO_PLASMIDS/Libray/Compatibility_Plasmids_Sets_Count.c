#include <MODEL.h>

int Compatibility_Plasmid_Sets_Count( Parameter_Table * Table, 
                                      int * Plasmid_Set_0, int N_0, int * Plasmid_Set_1, int N_1 )
{ 
  /* This function calculates the number of plasmids in the Donor (Plasmid_Set_0) that are fully
     compatible with all the plasmids in the recipient (Plasmid_Set_1)
  */       
          int k, l; 
          int n_Comp, n_Count;  
          int pd, pr;

          n_Count = 0; 
          for(k = 0; k < N_0; k++) {
            n_Comp = 0; 
            pd = Plasmid_Set_0[k];

            for(l = 0; l < N_1; l++) {
              pr = Plasmid_Set_1[l];
              if( Table->CPP[pd][pr] == 1.0 ) n_Comp++; 
            }

            if(n_Comp == N_1) n_Count++; 
          }

          return(n_Count);
}
