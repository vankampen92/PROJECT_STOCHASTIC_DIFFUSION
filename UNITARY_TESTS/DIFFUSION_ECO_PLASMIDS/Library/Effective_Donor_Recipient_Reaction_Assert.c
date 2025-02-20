#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

int Effective_Donor_Recipient_Reaction_Assert(Parameter_Table * Table, 
                                              int Donor_ID, int Donor_Sp, int * Plasmid_Donor_Potential_Set, int n_D, 
                                              int Recip_ID, int Recip_Sp, int k_R)
{ 
  /* During donor-recipient conjugation, a plasmid set from the donor is randomly selected from its 
     whole plasmid set. For the effective conjugation to occur, this assert should be passed. 
     This assert verifies if a conjugating plasmid set from the donor, this is, a given subset from 
     the whole donor plasmid profile, effectively changes the plasmid profile of the recipient. 
     This requires: 

     1. Permitted conjugation between (bacteria) strain types (conjugation condition)
     2. The full conjugating set (randomly selected from the donor) should be able to infect the recipient 
        bacterial type (full infection condition).
     3. The full conjugating set (randomly selected from the donor) should be totally compatible with all  
        plasmids in the recipient (full compatibility condition). 
     4. There is at least one plasmid in the conjugation plasmid set from the donor that is missing in the 
        recipient (weak complementarity condition). 
  */
    int N, k, n_bool; 
    int * Plasmid_Recipient_Set;  

    n_bool = 0; 
    if(Table->HBB[Recip_Sp][Donor_Sp] > 0.0) {   /* Conjugation is permitted */ 

      Plasmid_Recipient_Set = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int));  

      N = 0;
      for(k = 0; k<Table->No_of_PLASMIDS; k++)   
        if(Table->Strain_Profiles[Recip_Sp][k_R][k] == 1) { 
          /* Plasmid k is present in the RECIPIENT */
          Plasmid_Recipient_Set[N] = k; 
          N++;
        }
      Plasmid_Recipient_Set[Table->No_of_PLASMIDS] = N;

      N = 0;
      for(k = 0; k<n_D; k++) /* Full infection condition is required */
        if(Table->IBP[Recip_Sp][Plasmid_Donor_Potential_Set[k]] > 0.0) 
          /* Plasmid k-th plasmid potential donor set can infect potential recipient of bacterial type Recip_Sp */ 
          N++;
          
      if (N == n_D) {
        /* Full infection condition is fulfilled */
        n_bool = Compatibility_Plasmid_Sets_Assert( Table, 
                                                    Plasmid_Donor_Potential_Set, n_D, 
                                                    Plasmid_Recipient_Set, Plasmid_Recipient_Set[Table->No_of_PLASMIDS]);
        if( n_bool == 1) {
          /* Full compatibility condition is fulfilled */
          N = 0;
          for(k = 0; k<n_D; k++) 
             if(Table->Strain_Profiles[Recip_Sp][k_R][Plasmid_Donor_Potential_Set[k]] == 0)  
               /* Plasmid k-th from potential donor set is not present in potencial recipient strain */ 
               N++;
                
          if( N > 0 ) n_bool = 1;  /* Weak complementarity condition is fulfilled */
          else        n_bool = 0; 
        }                                        
      }      

      free(Plasmid_Recipient_Set);
    }

    return(n_bool);     
}
