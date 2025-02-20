#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

int Determining_actual_No_of_RESOURCES(Parameter_Table * Table)
{
  int N; /* Actual No of viable RESOURCES (actual number of different subpopulations)*/
  
  int n_bool;
  int No_of_POTENTIAL_PROFILES_per_STRAIN;
  int i, j, k; 

  No_of_POTENTIAL_PROFILES_per_STRAIN = (int)pow(2.0, Table->No_of_PLASMIDS);

  int ** Profile = (int **)calloc(No_of_POTENTIAL_PROFILES_per_STRAIN, sizeof(int *));
  for(i=0; i < No_of_POTENTIAL_PROFILES_per_STRAIN; i++) 
    Profile[i] = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int));

  /* No of actual viable profiles per strain. It can be different per each strain */ 
  int *  n = Table->n;  
  int *** Strain_Profiles = Table->Strain_Profiles; 

  Create_Binary_Combination( Profile, No_of_POTENTIAL_PROFILES_per_STRAIN, Table->No_of_PLASMIDS );

  for(i = 0; i < No_of_POTENTIAL_PROFILES_per_STRAIN; i++) {
  
    printf("%d:\t\t Binary Combination[%d] = [ ", i, i);
    for(j = 0; j<Table->No_of_PLASMIDS; j++)
      printf("%d ", Profile[i][j]);
    printf("]\n");
  }
  printf("\n");
  getchar();

  N = 0;  
  for(i = 0; i < Table->No_of_STRAINS; i++) {

    Table->n_0[i] = N;  /* Saving the Strain ID of a plasmid free strain */

    n[i] = 0;  
    for(j = 0; j<No_of_POTENTIAL_PROFILES_per_STRAIN; j++) {

      n_bool = Profile_Selfconsistency_Assert(i, Profile[j], Table);

      if(n_bool == 1) { 

        for(k = 0; k<Table->No_of_PLASMIDS + 1; k++) {
          if( k < Table->No_of_PLASMIDS )         
            Strain_Profiles[i][n[i]][k] = Profile[j][k];
          else
            Strain_Profiles[i][n[i]][k] = SumUP_Profile ( Profile[j], Table->No_of_PLASMIDS );
        }

        Table->StrainType_and_Profile[N][0] = i;     /* Strain Type (for Strain ID 'N')          */
        Table->StrainType_and_Profile[N][1] = n[i];  /* Profile No Specifier (for Strain ID 'N') */

        n[i]++;
        N++;                
      }      
    }
  }            
  Table->No_of_RESOURCES = N; 

  printf("Number of different strains and profiles: %d \n", N);
  N = 0; 
  for(i = 0; i < Table->No_of_STRAINS; i++) {

    assert( Table->n_0[i] == N);
    
    for(j = 0; j < n[i]; j++) {
      printf("Strain ID [%d]: Bacterial Type [%d]\t Profile: ", N, i );

      Print_Infection_Profile (Strain_Profiles[i][j], i, Table);

      printf("[ ");      
      for(k=0; k<Table->No_of_PLASMIDS; k++)
        printf("%d ", Strain_Profiles[i][j][k]);

      printf("]\n");

      N++;  
    }
    printf("\n");
  }  

  getchar(); 

  assert( N == Table->No_of_RESOURCES );

  for(i=0; i < No_of_POTENTIAL_PROFILES_per_STRAIN; i++) 
    free(Profile[i]);
  free(Profile);

  return(N); 
}
