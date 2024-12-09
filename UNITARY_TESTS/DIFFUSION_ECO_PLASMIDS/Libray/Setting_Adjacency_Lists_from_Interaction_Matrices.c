#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Setting_Adjacency_Lists_from_Interaction_Matrices (Parameter_Table *Table)
{
  int i, j, k, l, N, n_bool;
  int nc, nh, ni, nd;  
  int i_Focal, i_List, i_List_0;
  int k_Focal, k_List, k_List_0; 
  
  /* Adjacent Lists: */
  Table->No_of_CONJUGATION_EVENTS = 0;

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
    Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal);
    i_List_0 = Table->StrainType_and_Profile[i][0];
    k_List_0 = Table->StrainType_and_Profile[i][1];
    assert(i_Focal == i_List_0);
    assert(k_List_0 == k_Focal);

    assert( i_Focal >= 0 && i_Focal < Table->No_of_STRAINS);
    assert( k_Focal >= 0 && k_Focal < Table->n[i_Focal]);

    nc = 0; nh = 0; 
    for(j = 0; j<Table->No_of_RESOURCES; j++) {
      Calculate_Strain_and_Profile(Table, j, &i_List, &k_List);
      if( Table->ABB[i_Focal][i_List] > 0.0 ) {
        Table->Competition_List_Indeces[i][nc] = j;
        Table->Competition_Induced_Death[i][nc] = Table->ABB[i_Focal][i_List];
        nc++;
      }
      if( Table->HBB[i_Focal][i_List] > 0.0 )
        Table->Conjugation_List_Indeces[i][nh++] = j;      
    }
    Table->Competition_List_Indeces[i][Table->No_of_RESOURCES] = nc;     
    Table->Conjugation_List_Indeces[i][Table->No_of_RESOURCES] = nh;
    
    /* Recipiet list of the strain with the true index i (the focal DONOR) */
    N = 0;
    for(k = 0; k<Table->No_of_PLASMIDS; k++)
     if(Table->Strain_Profiles[i_Focal][k_Focal][k] == 1) 
      N++;
     
    if(N == 0) { 
      printf("Empty recipient list of focal donor (Strain ID: %d\t Strain Type: %d): ", i, i_Focal);
      Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Focal], Table->No_of_PLASMIDS); printf("\n"); 
      assert( k_Focal == 0 );
      Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES] = 0;
    }
    else {
      printf("List of RECIPIENTS of focal DONOR (Strain ID: %d\t Strain Type: %d): ", i, i_Focal);
      Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Focal], Table->No_of_PLASMIDS); printf("\n");

      ni = 0; 
      for(j = 0; j<Table->No_of_RESOURCES; j++) 
        if( j != i ) {
          Calculate_Strain_and_Profile(Table, j, &i_List, &k_List);  
          i_List_0 = Table->StrainType_and_Profile[j][0];           /* (i_List, k_List):   Potential RECIPIENT */
          k_List_0 = Table->StrainType_and_Profile[j][1];           /* (i_Focal, k_Focal): Focal DONOR )       */
          assert(i_List == i_List_0);
          assert(k_List_0 == k_List);

          n_bool = Potential_Donor_Recipient_Pair_Assert(Table, 
                                                         i, i_Focal, k_Focal,   /* (i_Focal, k_Focal): Focal DONOR )      */
                                                         j, i_List, k_List );   /* (i_List, k_List):  Potential RECIPIENT */
          if(n_bool == 1) {
            /* There is at least one plasmid in the donor potential set that is perfectibly 
               compatible with the whole plasmit set in the recipient, can infect that bacterial type 
               and can potentially change the plasmid profile of the recipieint 
            */
            Table->Recipient_List_Indeces[i][ni++] = j;   
          }
        }  
      Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES] = ni;
      Table->n_R[i]                                            = ni; /* No of Recipients of Strain ID 'i' */
      Table->No_of_CONJUGATION_EVENTS                         += ni;  
      printf("\n\n");
      // getchar();     
    }
    
    /* Donor list of a strain with the true index Strain ID 'i' (the focal RECIPIENT) */
    N = 0;
    for(k = 0; k<Table->No_of_PLASMIDS; k++)
     if(Table->Strain_Profiles[i_Focal][k_Focal][k] == 1) 
      N++;  

    if(N == Table->No_of_PLASMIDS) {
      printf("Empty DONOR list of FOCAL recipient (Strain ID: %d\t Strain Type: %d): ", i, i_Focal);
      Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Focal], Table->No_of_PLASMIDS); printf("\n");
      Table->Donor_List_Indeces[i][Table->No_of_RESOURCES] = 0;
    }
    else {
      printf("List of DONORS of a focal RECIPIENT (Strain ID: %d\t Strain Type: %d): ", i, i_Focal);
      Print_Strain_Profile(Table->Strain_Profiles[i_Focal][k_Focal], Table->No_of_PLASMIDS); printf("\n");

      nd = 0; 
      for(j = 0; j<Table->No_of_RESOURCES; j++) 
        if( j != i ) { 
          Calculate_Strain_and_Profile(Table, j, &i_List, &k_List);   /* Potential Donor */
          i_List_0 = Table->StrainType_and_Profile[j][0];
          k_List_0 = Table->StrainType_and_Profile[j][1];
          assert(i_List   == i_List_0);
          assert(k_List_0 == k_List); 

          n_bool = Potential_Donor_Recipient_Pair_Assert(Table,   
                                                         j, i_List, k_List ,      /* (i_List, k_List):   Potential DONOR */
                                                         i, i_Focal, k_Focal );   /* (i_Focal, k_Focal): Focal RECIPIENT */
          if(n_bool == 1) {
            /* There is at least one plasmid in the donor potential set that is perfectibly 
               compatible with the whole plasmit set in the recipient, can infect that bacterial type, and can 
               change the recipient's profile 
            */
            Table->Donor_List_Indeces[i][nd++] = j;   
          }
        }       
        
      Table->Donor_List_Indeces[i][Table->No_of_RESOURCES] = nd;
      printf("\n\n");
      // getchar();
    } 
  }

  /* Plasmid Compatitibility List */
  for(i = 0; i<Table->No_of_PLASMIDS; i++) {
    N = 0; 
    for(j = 0; j<Table->No_of_PLASMIDS; j++) {
      if(Table->CPP[i][j] > 0.0)
        Table->Plasmid_Compatibility_Indeces[i][N++] = j;
    }
    Table->Plasmid_Compatibility_Indeces[i][Table->No_of_PLASMIDS] = N;
  }  
}
