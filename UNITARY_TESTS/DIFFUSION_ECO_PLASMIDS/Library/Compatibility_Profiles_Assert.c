#include <MODEL.h>

int Compatibility_Profiles_Assert (Parameter_Table * Table, int * Profile_0, int * Profile_1, int N) 
{ 
  int k; 
  int n_bool; 
  int N_0, N_1;     

  int * Plasmid_Set_0 = (int *)calloc(N, sizeof(int));
  int * Plasmid_Set_1 = (int *)calloc(N, sizeof(int));

          N_0 = 0; 
          for(k = 0; k<N; k++) 
            if( Profile_0[k] == 1)
              Plasmid_Set_0[N_0++] = k; 
           
          N_1 = 0; 
          for(k = 0; k<N; k++) 
            if( Profile_1[k] == 1)
              Plasmid_Set_1[N_1++] = k;

          n_bool = Compatibility_Plasmid_Sets_Assert( Table, Plasmid_Set_0, N_0, Plasmid_Set_1, N_1);

  free(Plasmid_Set_0);
  free(Plasmid_Set_1);

  return(n_bool);    
} 
