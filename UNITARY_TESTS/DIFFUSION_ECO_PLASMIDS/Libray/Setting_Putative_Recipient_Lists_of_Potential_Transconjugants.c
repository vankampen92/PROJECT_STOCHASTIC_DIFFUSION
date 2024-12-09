#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Setting_Putatitive_Recipient_Lists_of_Potential_Trasconjugants (Parameter_Table * Table)
{
  int L, N, i, k, l, i_Strain_ID, k_Strain_ID, ni; 
  int i_Focal, k_Focal;  
  int * Strain_Profile; 
  int * Pl_T; 
  int ** Profile; 

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
    Pl_T = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int)); /* TRANSCONJUGANT  */

    Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal);

    /* Plasmid set of strain with the true Strain ID i (the focal TRANSCONJUGANT) */
    N = 0;
    for(k = 0; k<Table->No_of_PLASMIDS; k++)
     if(Table->Strain_Profiles[i_Focal][k_Focal][k] == 1) {
      Pl_T[N] = k;
      N++;
     }
    Pl_T[Table->No_of_PLASMIDS] = N;  /* The Cardinal of this Set */  
    assert( N == Table->Strain_Profiles[i_Focal][k_Focal][Table->No_of_PLASMIDS] );

    if( N > 0 ) {
      printf("Focal Transconjugant: (Strain ID: %d) (Strain TYPE: %d): ", i, i_Focal);
      Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Focal], Table->No_of_PLASMIDS); printf("\n");

      L = (int)pow( 2.0, (double)N );

      Profile = (int **)calloc(L, sizeof(int *));
      for(k=0; k < L; k++) 
        Profile[k] = (int *)calloc(N, sizeof(int));

      Create_Binary_Combination( Profile, L, N );
    
      ni = 0; 
      for(k = 0; k < L-1; k++) {    
        Strain_Profile = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int));
      
        for(l=0; l<Pl_T[Table->No_of_PLASMIDS]; l++) 
          Strain_Profile[Pl_T[l]] = Profile[k][l];

        Strain_Profile[Table->No_of_PLASMIDS] = SumUP_Profile (Strain_Profile, Table->No_of_PLASMIDS);       
        Print_Strain_Profile(Strain_Profile, Table->No_of_PLASMIDS); printf("\n");

        i_Strain_ID = Calculate_Strain_ID_from_Profile_and_Strain_Type( Table, 
                                                                       i_Focal, Strain_Profile );
        if( i_Strain_ID == Table->No_of_RESOURCES ) {
          Printing_Strains_Profiles(Table);

          printf("Impossible Strain ID (%d)!!!\n", i_Strain_ID);
          printf("\n ------------------------------------------------------------------------------ \n");
          printf("No of RESOURCES (No of STRAIN IDs): %d\n", Table->No_of_RESOURCES);
          printf("No of STRAIN TYPES: %d\n", Table->No_of_STRAINS);
          printf("No of PLASMIDS: %d\n", Table->No_of_PLASMIDS);
          printf("No of TOTAL EVENTS: %d\n", Table->No_of_STRAINS * 6 + Table->No_of_CONJUGATION_EVENTS);
          printf("No of Conjugation Pairs (Donor, Recipient): %d\n", Table->No_of_CONJUGATION_EVENTS);
          
          exit(0); 
        } 

        assert(i_Strain_ID >= Table->n_0[i_Focal]);
        assert(i_Strain_ID < Table->n_0[i_Focal] + Table->n[i_Focal]);                                                               

        Table->Putative_Recipient_List_Indeces[i][ni++] = i_Strain_ID;

        k_Strain_ID = Table->StrainType_and_Profile[i_Strain_ID][1];   

        if (Table->StrainType_and_Profile[i_Strain_ID][0] != i_Focal )
          printf("Warning: i_Strain_ID = %d\t i_Focal = %d\n", 
                  Table->StrainType_and_Profile[i_Strain_ID][0], i_Focal);

        assert(Table->StrainType_and_Profile[i_Strain_ID][0] == i_Focal);

        // Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Strain_ID], Table->No_of_PLASMIDS); printf("\n");
        // printf("\n");

        free(Strain_Profile);                                                        
      }
      Table->Putative_Recipient_List_Indeces[i][Table->No_of_PROFILES] = ni;

      for(k=0; k < L; k++) 
        free(Profile[k]);
      free(Profile);
    }
    else {
      printf("Focal Potential Transconjugant: (Strain ID: %d) (Strain TYPE: %d): ", i, i_Focal);
      Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Focal], Table->No_of_PLASMIDS); printf("\n");
      printf(" Plasmid-free profiles can not be transconjugants!!!\n");
      printf("They should have a null putative recipient list.\n");
      Table->Putative_Recipient_List_Indeces[i][Table->No_of_PROFILES] = 0;
    }

    // getchar();  
    free(Pl_T);
  }  
}
