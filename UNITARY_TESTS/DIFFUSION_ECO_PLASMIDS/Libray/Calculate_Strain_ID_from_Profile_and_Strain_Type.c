#include <MODEL.h>

int Calculate_Strain_ID_from_Profile_and_Strain_Type(Parameter_Table * Table, 
                                                     int Sp_Strain, int * Profile )
{
  /* Input:
     . Table
     . Sp_Strain, Strain Type: Sp_Strain >= 0 && Sp_Strain < No_of_STRAINS
     . Profile, Potential profile to search for. 

     Output:
     . Strain ID corresponding to input Strain Type 'Sp_Strain' and 'Profile'  
  */
  int i, N, N_0, k, n; 
  int Sp_Strain_ID;
  int n_Bool_Full; 

  N_0 = Table->n_0[(Sp_Strain)];

  N = 0; n = 0; n_Bool_Full = 0;  
  while( N < Table->n[Sp_Strain] && n < Table->No_of_PLASMIDS ) {

    n = 0;
    for( k = 0; k < Table->No_of_PLASMIDS; k++)
      if(Table->Strain_Profiles[Sp_Strain][N][k] == Profile[k])
        n++; 

    if( n == Table->No_of_PLASMIDS ) n_Bool_Full = 1; 

    N++;    
  }

  if( n_Bool_Full == 0 ) {
    Sp_Strain_ID = Table->No_of_RESOURCES;       /* Impossible Strain ID: No profile found */
  }
  else {
    Sp_Strain_ID =  N_0 + N - 1;                 /* Strain_ID associated to the input Profile[] and input Sp_Strain Type */
  }

  return(Sp_Strain_ID);
}
