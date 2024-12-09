#include <MODEL.h>

int Compatibility_Plasmid_Sets_Assert( Parameter_Table * Table, 
                                       int * Plasmid_Set_0, int N_0, int * Plasmid_Set_1, int N_1 )
{ 
          int k, l; 
          int n_bool; 
          int n_Comp;  
          int pd, pr;

          n_Comp = 0; 
          for(k = 0; k < N_0; k++) {
            pd = Plasmid_Set_0[k];
            for(l = 0; l < N_1; l++) {
              pr = Plasmid_Set_1[l];
              if( Table->CPP[pd][pr] == 1.0 ) n_Comp++; 
            }   
          }
          if(n_Comp == N_0*N_1) n_bool = 1;
          else                  n_bool = 0; 

          return(n_bool);
}
