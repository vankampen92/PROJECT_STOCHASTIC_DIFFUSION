#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

int Potential_Donor_Recipient_Pair_Assert(Parameter_Table * Table, 
                                          int Donor_ID, int Donor_Sp, int k_D,
                                          int Recip_ID, int Recip_Sp, int k_R)
{
  /* This verifies if the plasmid set from the donor can potentially change the plasmid profile 
     of the recipient. This requires:

     1. Permitted conjugation between bacterial types (conjugation condition)
     2. At least one plasmid from the donor set should be able to infect the recipient bacterial type
        (weak infetion condition) 
     3. From those plasmids that can infect, there is at least one plasmid that is missing in the recipient 
        (weak complementarity condition).
     4. From those plasmids that can infect and are missing in the recipient, at least one should be totally 
        compatible with all plasmids in the recipient (weak compatibility condition). 
  */
    int N, k, n_bool, n_Count; 
    int * Plasmid_Recipient_Set;
    int * Plasmid_Donor_Potential_Set;    

    n_bool = 0; 
    if(Table->HBB[Recip_Sp][Donor_Sp] > 0.0) {      /* Conjugation is permitted */ 
      
      Plasmid_Recipient_Set        = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int));   
      Plasmid_Donor_Potential_Set  = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int));   

      N = 0;
      for(k = 0; k<Table->No_of_PLASMIDS; k++)   
        if(Table->Strain_Profiles[Recip_Sp][k_R][k] == 1) { 
          /* Plasmid k is present in the RECIPIENT */
          Plasmid_Recipient_Set[N] = k; 
          N++;
        }
      Plasmid_Recipient_Set[Table->No_of_PLASMIDS] = N;

      N = 0;
      for(k = 0; k<Table->No_of_PLASMIDS; k++) /* Infection is permitted */
        if(Table->Strain_Profiles[Donor_Sp][k_D][k] == 1 && Table->IBP[Recip_Sp][k] > 0.0) {
          /* Plasmid k from focal donor Strain ID 'i' can infect potential recipient Strain ID 'j' */ 
          if(Table->Strain_Profiles[Recip_Sp][k_R][k] == 0) { 
            /* Plasmid k from focal donor Strain ID 'i' is not present in potencial recipient Strain ID 'j' */
            Plasmid_Donor_Potential_Set[N] = k; 
            N++;
          }    
        }
      Plasmid_Donor_Potential_Set[Table->No_of_PLASMIDS] = N;

      if( N > 0) {
        n_Count = Compatibility_Plasmid_Sets_Count(Table,  
                                                 Plasmid_Donor_Potential_Set, Plasmid_Donor_Potential_Set[Table->No_of_PLASMIDS], 
                                                 Plasmid_Recipient_Set, Plasmid_Recipient_Set[Table->No_of_PLASMIDS]);
        if ( n_Count > 0 ) n_bool = 1;

        if(n_bool == 1) {
          /* There is at least one plasmid in the donor potential set that is perfectibly 
            compatible with the whole plasmit set in the recipient, can infect that bacterial type, and can 
            change the recipient's profile 
          */
          printf("Donor: "); 
          printf(" Y\t"); Print_Infection_Profile(Table->Strain_Profiles[Donor_Sp][k_D], Recip_Sp, Table); printf("\t"); 
          printf("Recipient: "); 
          Print_Strain_Profile(Table->Strain_Profiles[Recip_Sp][k_R], Table->No_of_PLASMIDS); printf("\n\n");
        }
        else {
          printf("Donor: "); 
          printf(" N\t"); Print_Infection_Profile(Table->Strain_Profiles[Donor_Sp][k_D], Recip_Sp, Table); printf("\t"); 
          printf("Recipient: "); 
          Print_Strain_Profile(Table->Strain_Profiles[Recip_Sp][k_R], Table->No_of_PLASMIDS); printf("\n\n"); 
        } 
      }

      free(Plasmid_Recipient_Set);
      free(Plasmid_Donor_Potential_Set);              
    }

  return(n_bool);
}
