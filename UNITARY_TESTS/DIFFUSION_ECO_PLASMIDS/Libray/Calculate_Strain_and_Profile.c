#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Calculate_Strain_and_Profile(Parameter_Table * Table, int n, int * i_Strain, int * k_Profile)
{
  /* Input: 
     . n, label of subpopulation or local variable index (Strain ID or true index of the local variable)
 
     Output:
     . i_Strain, label of the strain
     . k_Profile, label of the profile (within the i_Strain)
  */

  int i;
  int N; 
  int Sp_Strain, Sp_Profile;   
   
  N = Table->n[0]; 
  i = 0; 
  while ( (n - N) >= 0 )
  {
    N += Table->n[i+1];
    i++;   
  }
  * i_Strain = i; 
  
  N = 0; 
  for(i=0; i < (* i_Strain); i++)
  { 
    N += Table->n[i]; 
  }
  * k_Profile = n - N; 
  
  // * i_Strain  = Table->StrainType_and_Profile[n][0];
  // * k_Profile = Table->StrainType_and_Profile[n][1];

  Sp_Strain  = Table->StrainType_and_Profile[n][0];
  Sp_Profile = Table->StrainType_and_Profile[n][1];

  assert( (* i_Strain) == Sp_Strain && (* k_Profile) == Sp_Profile );

}
