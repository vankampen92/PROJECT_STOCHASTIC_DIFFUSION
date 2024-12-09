#include <MODEL.h>

void Setting_Reactive_Recipient_Donor_Pairs_and_Rates(Parameter_Table * Table) 
{
  ///   A potential conjugation process may fail if the selected plasmids to transmit from the donor fail 
  ///   to pass the full infection conditiion or the full compatibility condition or, at least, the (weak) 
  ///   complementarity condition.
  int N;
  int k, l, n; 
  int n_Effective_Pair; 
  int * Donor_List_Subset;
  double * Donor_Rate_Subset; 
  double Rate;
  int Strain_ID_R, Strain_ID_D;
  int * Profile;   

    for(k = 0; k < Table->No_of_RESOURCES; k++) {
        for(l = 0; l < Table->Putative_Recipient_List_Indeces[k][Table->No_of_PROFILES]; l++) { 

          Strain_ID_R = Table->Putative_Recipient_List_Indeces[k][l]; 

          Donor_List_Subset = (int *)calloc(Table->Donor_List_Indeces[Strain_ID_R][Table->No_of_RESOURCES], 
                                            sizeof(int));
          Donor_Rate_Subset = (double *)calloc(Table->Donor_List_Indeces[Strain_ID_R][Table->No_of_RESOURCES], 
                                               sizeof(double));                                  
          N = 0; 
          for(n = 0; n < Table->Donor_List_Indeces[Strain_ID_R][Table->No_of_RESOURCES]; n++) {
            Strain_ID_D = Table->Donor_List_Indeces[Strain_ID_R][n];  

            /* . k:           Transconjugant Strain ID */
            /* . Strain_ID_R: Recipient Strain ID      */
            /* . Strain_ID_D: Donor Strain ID          */
            n_Effective_Pair = Recipient_Donor_Transconjugant_Rate( k, Strain_ID_R, Strain_ID_D, 
                                                                    &Rate, Table ); 
            if(n_Effective_Pair == 1) {
              Donor_List_Subset[N] = Strain_ID_D; 
              Donor_Rate_Subset[N] = Rate;
              N++; 
            }
            else{
              printf(" Effective Conjugation would never produce the required transconjugant:\n");
              
              Profile = Table->Strain_Profiles[Table->StrainType_and_Profile[k][0]][Table->StrainType_and_Profile[k][1]];
              Print_Strain_Profile(Profile, Table->No_of_PLASMIDS); 
              printf("\t Transconjugant Profile");; printf("\n");
              
              Profile = Table->Strain_Profiles[Table->StrainType_and_Profile[Strain_ID_D][0]][Table->StrainType_and_Profile[Strain_ID_D][1]];
              Printf_Infection_Profile(Profile, Table->StrainType_and_Profile[Strain_ID_R][0], Table); 
              printf("\t\t Donor Profile (when infecting the recipient)"); printf("\n");
              
              Profile = Table->Strain_Profiles[Table->StrainType_and_Profile[Strain_ID_R][0]][Table->StrainType_and_Profile[Strain_ID_R][1]]; 
              Print_Strain_Profile(Profile, Table->No_of_PLASMIDS);
              printf("\t Recipient Profile."); printf("\n\n"); 
              // getchar(); 
            }            
          }

          Table->DoRe[k][l] = (Donor_Recipient_Pair **)calloc(N, sizeof(Donor_Recipient_Pair *));              
          for( n = 0; n < N; n++) {
            Table->DoRe[k][l][n] = (Donor_Recipient_Pair *)calloc(1, sizeof(Donor_Recipient_Pair));
            Table->DoRe[k][l][n]->ne = Donor_List_Subset[n];
            Table->DoRe[k][l][n]->ra = Donor_Rate_Subset[n];
          
            Table->DoRe[k][l][n]->N  = N; 
          }  

          // printf(" Number of effective donors of the %d-th putative recipient of strain ID %d:\t %d (out of %d)\n", l, k,  
          //          Table->DoRe[k][l][0]->N, 
          //          Table->Donor_List_Indeces[Strain_ID_R][Table->No_of_RESOURCES]); 

          free(Donor_List_Subset);
          free(Donor_Rate_Subset);                                                    
        }
    }
}
