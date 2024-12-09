#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void show_DoubleMatrix(double **M, int Nx, int Ny)
{
  int i,j, n_sys;

  /*n_sys = system("clear");*/
  printf("Printing matrix...\n");

  printf("* :\t");
  for(i=0; i<Ny; i++) printf("%4d  ",i);
  printf("\n");
  printf("   \t");
  for(i=0; i<Ny; i++) printf("..... ");
  printf("\n");

  for(i=0; i<Nx; i++){
    printf("%d :\t",i);
    for(j=0; j<Ny; j++) printf("%4.3g  ",M[i][j]);
    printf("\n");
  }
  printf("\n\n");
  /* getchar(); */
}

void Setting_Interaction_Matrices (Parameter_Table * Table) 
{
    int i, j;

    // assert(Table->TYPE_of_MODEL == 21);  /* 21: MODEL=DIFFUSION_ECO_PLASMIDS */

    for (i=0; i<Table->No_of_STRAINS; i++) {
      
      for (j=i; j<Table->No_of_STRAINS; j++) {
        // ABB, competition matrix: symmetric competition

        if(j == i) 
          Table->ABB[i][j] = 0.0;

        else if(gsl_rng_uniform_pos(r) < Table->p_2){
          
          Table->ABB[i][j] = Table->Delta_C_0;
          Table->ABB[j][i] = Table->Delta_C_0;
        }
        
        else {
          Table->ABB[i][j] = 0.0;
          Table->ABB[j][i] = 0.0;
        } 
      }

      for (j=i; j<Table->No_of_STRAINS; j++) {
        // HBB Conjugation Matrix  
        
        if(j == i) 
          Table->HBB[i][j] = 1.0;      /* Self-conjugation between two cells of the 
                                          same bacterial type is always possible.
                                          Efective conjungation will depend on their
                                          respective plasmid profiles.  
                                          However, conjugation between two individual 
                                          cells with the same plasmid profile never 
                                          yields a transconjungant different from the 
                                          recipient. It is an event that will not 
                                          change the configuration of the system. 
                                          It will not be considered. 
                                       */
        else if (gsl_rng_uniform_pos(r) < Table->p_2){
           Table->HBB[i][j] = 1.0;
           Table->HBB[j][i] = 1.0;
        }

        else {
          Table->HBB[i][j] = 0.0;
          Table->HBB[j][i] = 0.0; 
        }
      }
         
      for (j=0; j<Table->No_of_PLASMIDS; j++) {
        // IBP Bacteria-Plasmid Interaction Matrix  (Bipartite Network)
        if(gsl_rng_uniform_pos(r) < Table->p_2)
          Table->IBP[i][j] = 1.0;
        else
          Table->IBP[i][j] = 0.0;
      }
    }

    // CPP, Plasmid-Plasmid Compatibility Matrix
    for (i=0; i<Table->No_of_PLASMIDS; i++) {

      for (j=i+1; j<Table->No_of_PLASMIDS; j++) {
        if(gsl_rng_uniform_pos(r) < Table->p_2) {
          Table->CPP[i][j] = 1.0; 
          Table->CPP[j][i] = 1.0;
        }
        else {
          Table->CPP[i][j] = 0.0;
          Table->CPP[j][i] = 0.0;
        }
      }

      Table->CPP[i][i] = 1.0; /* Reinforicing plasmid self-compatibility */
    } 
}

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

      Printf_Infection_Profile (Strain_Profiles[i][j], i, Table);

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

int SumUP_Profile ( int * Profile, int N )
{
  int i, S;

  S = 0; 
  for( i=0; i<N; i++)
    S += Profile[i];

  return(S);  
}

int Infection_Condition_Assert(int i_Sp, int * Profile, Parameter_Table * Table)
{
  int k; 
  int n_bool; 
  int n_count; 
  int S; 

  S = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++) 
      if( Profile[k] == 1) S++;

    /* Counting infections */
  n_count = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++) {
      if( Profile[k] == 1) {  
            /* Infection Constraint */ 
            if(Table->IBP[i_Sp][k] == 1.0) n_count++; 
      }
    }
    
  if(n_count == S) n_bool = 1;

  return(n_bool); 
}

int Profile_Selfconsistency_Assert(int i_Strain_Sp, int * Profile, Parameter_Table * Table)
{
  int assert_true; 
  int assert_infectivity; 
  int assert_compatibility; 
  int S, k; 

    S = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++) 
      if( Profile[k] == 1) S++;

    if (S == 0) { /* Plasmid-free Profile */  
          assert_true = 1; 
          return(assert_true);
    }
    else {
          assert_infectivity = Infection_Condition_Assert(i_Strain_Sp, Profile, Table); 

          assert_compatibility = Compatibility_Profiles_Assert (Table, Profile, Profile, 
                                                                Table->No_of_PLASMIDS); 

          if( assert_infectivity == 1 && assert_compatibility == 1) assert_true = 1;
          else                                                      assert_true = 0; 

          return(assert_true);
    }
}

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
          printf(" Y\t"); Printf_Infection_Profile(Table->Strain_Profiles[Donor_Sp][k_D], Recip_Sp, Table); printf("\t"); 
          printf("Recipient: "); 
          Print_Strain_Profile(Table->Strain_Profiles[Recip_Sp][k_R], Table->No_of_PLASMIDS); printf("\n\n");
        }
        else {
          printf("Donor: "); 
          printf(" N\t"); Printf_Infection_Profile(Table->Strain_Profiles[Donor_Sp][k_D], Recip_Sp, Table); printf("\t"); 
          printf("Recipient: "); 
          Print_Strain_Profile(Table->Strain_Profiles[Recip_Sp][k_R], Table->No_of_PLASMIDS); printf("\n\n"); 
        } 
      }

      free(Plasmid_Recipient_Set);
      free(Plasmid_Donor_Potential_Set);              
    }

  return(n_bool);
}

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

void Printing_Putative_Recipient_Lists(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 
  int i_List, k_List; 

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
   Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 

    printf("\n\n");
    printf("Putative Recipient List of Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", i, i_Focal, k_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]:\n");   

    for(j=0; j<Table->Putative_Recipient_List_Indeces[i][Table->No_of_PROFILES]; j++) {

      Calculate_Strain_and_Profile(Table, Table->Putative_Recipient_List_Indeces[i][j], &i_List, &k_List);   

      printf("Strain ID %d:\t Bacterial Type (%d) and Profile (%d): [ ", 
        Table->Putative_Recipient_List_Indeces[i][j], i_List, k_List);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_List][k_List][k]);
      printf("]\n"); 
    }

    printf("\n");
    printf("Putatative Recipient List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Putative_Recipient_List_Indeces[i][Table->No_of_PROFILES]; j++)
      printf(" %d ", Table->Putative_Recipient_List_Indeces[i][j]);
    printf("}");

    printf("\n\n");
  }
}

void Printing_Strains_Profiles(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 

  for(i=0; i<Table->No_of_STRAINS; i++)
    printf("Strain Type: %d\t No of Profiles of this strain type: %d\n", i, Table->n[i]);
  printf("\n"); printf("\n");

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
    Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 
   
    printf("Strain ID: %d\t Strain Type: %d\t Strain Profile [ ", i, i_Focal);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]");

    printf("\n");  
  }
}

void Printing_Strains_Profiles_and_Lists(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 
  int i_List, k_List; 

  for(i=0; i<Table->No_of_STRAINS; i++)
    printf("Strain Type: %d\t No of Profiles of this strain type: %d\n", i, Table->n[i]);
  printf("\n");

  /* Adjacent Lists: */
  for(i = 0; i<Table->No_of_RESOURCES; i++) {
   Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 
   
   printf("Strain ID: %d\t Strain Type: %d\t Strain Profile [ ", i, i_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]");
  
    printf("\n");
    printf("Competition List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Competition_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Competition_List_Indeces[i][j]);
    printf("}"); 

    printf("\n");
    printf("Conjugation List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Conjugation_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Conjugation_List_Indeces[i][j]);
    printf("}");

    printf("\n\n");
    printf("Recipient List of Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", i, i_Focal, k_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]:\n");   
    for(j=0; j<Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES]; j++) {

      Calculate_Strain_and_Profile(Table, Table->Recipient_List_Indeces[i][j], &i_List, &k_List);      
      printf("Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", Table->Recipient_List_Indeces[i][j], i_List, k_List);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_List][k_List][k]);
      printf("]\n"); 
    }
    printf("\n");
    printf("Recipient List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Recipient_List_Indeces[i][j]);
    printf("}");

    printf("\n\n");

    printf("Donor List of Strain ID %d:\t Bacterial Type (%d) and Profile (%d): [ ", i, i_Focal, k_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]:\n");  
    for(j=0; j<Table->Donor_List_Indeces[i][Table->No_of_RESOURCES]; j++) {
      
      Calculate_Strain_and_Profile(Table, Table->Donor_List_Indeces[i][j], &i_List, &k_List);      
      printf("Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", Table->Donor_List_Indeces[i][j], i_List, k_List);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_List][k_List][k]);
      printf("]\n");      
    }
    printf("\n");
    printf("Donor List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Donor_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Donor_List_Indeces[i][j]);
    printf("}");
    printf("\n\n");
  }
  
  printf("\n");
  printf("Plasmid Compatibility List:\n");
  for(i=0; i<Table->No_of_PLASMIDS; i++) {
    printf("Plasmid ID: %d\t Compatibility list: [ ", i);
    for(j=0; j<Table->Plasmid_Compatibility_Indeces[i][Table->No_of_PLASMIDS]; j++)
      printf("%d ", Table->Plasmid_Compatibility_Indeces[i][j]);
    printf("]\n");
  }        
}

void Print_Strain_Profile(int * Profile, int No_of_PLASMIDS) 
{
    /* This input Profile should have been created to store
       (No_of_PLASMIDS + 1) integers 
    */
      int k;
      printf("[ ");
      for(k=0; k<No_of_PLASMIDS; k++) 
        printf("%d ", Profile[k]);
      printf("] ");
      printf("[[ %d ]]", Profile[No_of_PLASMIDS]);   
}

void Printf_Infection_Profile(int * Profile, int i_List, Parameter_Table * Table)
{
  int k;
      printf("[ ");
      
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        if(Profile[k] == 1 ) {
          if( Table->IBP[i_List][k] == 1.0) 
            printf("Y ");
          else
            printf("N ");
        }
        else
          printf("%d ", Profile[k]);

      printf("] ");  
}

void Create_Binary_Combination( int ** Binary_Combination, int N, int LENGTH )
{
  /* 
     This code creates an ordered list of binary numbers:
     from 1 to 2^{MAX_LENGTH}. Binary numbers are stored 
     in a int ** Binary_Combination array:
        
        0:                  Binary_Combination[0] = {00000000}
        1:                  Binary_Combination[1] = {00000001}
        2:                  Binary_Combination[2] = {00000010}
        .
.       .
        .
        2^{MAX_LENGTH}-1    Binary_Combination[2^{MAX_LENGTH}-1]

     This function creates exhaustively all 1/0 strings of a given 
     length LENGTH.      

     The numbers of strings to create should be given as the N input
     parameter. Enough space should have been reserved (in the parent function) in 
     Binary_Combination[][] to store all of them      
  */

    int * number  = (int *)calloc(LENGTH, sizeof( int ) );

    int_buffer_rec(Binary_Combination, N, number, LENGTH, LENGTH);

    free(number);
}

void int_buffer_rec(int ** Number_List, int N,
                    int * number, int n, int length)
{
    int i;
    static int m = 0;

    if(n > 0) {
        number[length - n] = 0;
        int_buffer_rec(Number_List, N, number, n - 1, length);
        number[length - n] = 1;
        int_buffer_rec(Number_List, N, number, n - 1, length);
    }
    else {
        for(i = 0; i < length; ++i) {
            // Rprintf("%u", number[i]);
            Number_List[m][i]=number[i];
        }
        // Rprintf("\n");
        m++;
    }

    if ( m == N ) m = 0;
}

void GSL_Init_Random_Seed(const gsl_rng * r)
{
  /* This function seed the GSL Random Generator r
     with a seed which is different for each initialization
     according to current computer time 
  */
        unsigned long int     seed;
        time_t  nowtime;
        struct  tm *preztime;

        time(&nowtime);
        preztime = localtime(&nowtime);
        seed = (int)((preztime->tm_sec+1)*(preztime->tm_min+1)*
                (preztime->tm_hour+1)*(preztime->tm_year)*(preztime->tm_year));
        if(seed%2==0) seed++;

        printf(" Random Number Seed: %lu\n", seed);

        gsl_rng_set(r, seed);
}

void GSL_Init_Random_Seed_from_File(const gsl_rng * r)
{
        unsigned long int     seed;

        //seed = 100;       
        /* A script have changed the environmetal variable
           GSL_RNG_SEED before the execution of the program
           starts. This value has been set to gsl_rng_default_seed
           in gsl_random_number_Setup.c. This setup is 
           always done at the start of code execution.
        */

        printf ("GSL_RNG_SEED = %lu\n", gsl_rng_default_seed);

        gsl_rng_set(r, gsl_rng_default_seed);

        printf ("first value = %lu\n", gsl_rng_get (r));
}

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

int Compatibility_Plasmid_Sets_Count( Parameter_Table * Table, 
                                      int * Plasmid_Set_0, int N_0, int * Plasmid_Set_1, int N_1 )
{ 
  /* This function calculates the number of plasmids in the Donor (Plasmid_Set_0) that are fully
     compatible with all the plasmids in the recipient (Plasmid_Set_1)
  */       
          int k, l; 
          int n_Comp, n_Count;  
          int pd, pr;

          n_Count = 0; 
          for(k = 0; k < N_0; k++) {
            n_Comp = 0; 
            pd = Plasmid_Set_0[k];

            for(l = 0; l < N_1; l++) {
              pr = Plasmid_Set_1[l];
              if( Table->CPP[pd][pr] == 1.0 ) n_Comp++; 
            }

            if(n_Comp == N_1) n_Count++; 
          }

          return(n_Count);
}

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

void Transconjugation_Gain_and_Loss_Total_Rates( const double * Y, double * Gain, double * Loss, 
                                                 Parameter_Table * Table )
{
  /* Input: 
        . Table
        . Y, state vector 
     Output
        . Gain[k]: production rate of transconjugants of Strain ID 'k'
        . Loss[k]: loss rate of recipients as they turn into transconjugants. 
  */
  int k, l, n, N;
  double G;
  int Strain_ID_R, Strain_ID_D;     
  double K_R = (double)Table->K_R; 
  
  double * L = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));   
  
  for(k = 0; k < Table->No_of_RESOURCES; k++) { 

    G = 0.0; 
    for(l = 0; l < Table->Putative_Recipient_List_Indeces[k][Table->No_of_PROFILES]; l++) {
        Strain_ID_R = Table->Putative_Recipient_List_Indeces[k][l];
        
        N = Table->DoRe[k][l][0]->N; 
        for(n = 0; n < N; n++) {
           
          Strain_ID_D     = Table->DoRe[k][l][n]->ne;    
          
          G              += Table->Eta_RP[Strain_ID_D] * Table->DoRe[k][l][n]->ra * Y[Strain_ID_R]/K_R * Y[Strain_ID_D];
          L[Strain_ID_R] += Table->Eta_RP[Strain_ID_D] * Table->DoRe[k][l][n]->ra * Y[Strain_ID_R]/K_R * Y[Strain_ID_D];
        }   
    }
    Gain[k] = G;
  }

  for(k = 0; k < Table->No_of_RESOURCES; k++) 
    Loss[k] = L[k];   
  
  free(L); 
}

int Recipient_Donor_Transconjugant_Rate( int Trans_ID, int Strain_ID_R, int Strain_ID_D, 
                                         double * Rate, Parameter_Table * Table )
{
  ///   This function calculates the rate of formation of the transconjugant, 'Trans_ID', for a given  
  ///   encounter between a Strain Type and Profile of the recipeint (specified by a given Strain_ID_R) 
  ///   and a donor (Strain_ID_D) as result of an efficient conjugation process.
  
  ///   Notice that this rate, (* Rate), is a probability, which will multiply the ecounter rate to build the 
  ///   effecitve rate of formation of the transconjugant as a result of the effective conjugation of a given 
  ///   reactive donor-recipient pair.  

  int n_bool; 
  int k, l, n, n_0, d, n_Difference_TR, n_Difference_DR, n_Invasion;
  int Strain_Sp_T, Strain_Sp_R, Strain_Sp_D; 
  int * Profile_T; 
  int * Profile_R; 
  int * Profile_D; 

  Profile_D = Table->Strain_Profiles[Table->StrainType_and_Profile[Strain_ID_D][0]][Table->StrainType_and_Profile[Strain_ID_D][1]];
  Profile_R = Table->Strain_Profiles[Table->StrainType_and_Profile[Strain_ID_R][0]][Table->StrainType_and_Profile[Strain_ID_R][1]];
  Profile_T = Table->Strain_Profiles[Table->StrainType_and_Profile[Trans_ID][0]][Table->StrainType_and_Profile[Trans_ID][1]];

  int * Pl_R = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int));  
  int * Pl_D = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int));
  int * Pl_Difference_TR = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int));  /* Set of Plasmids from the Transconjugant that
                                                                                  are NOT in the Recipient */
  int * Pl_Difference_DR = (int *)calloc(Table->No_of_PLASMIDS, sizeof(int));  /* Set of Plasmids from the Donor that
                                                                                  are NOT in the Recipient */
  n = 0; 
  for(k = 0; k<Table->No_of_PLASMIDS; k++)
    if (Profile_R[k] == 1)
      Pl_R[n++] = k; 
  n = 0; 
  for(k = 0; k<Table->No_of_PLASMIDS; k++)
    if (Profile_D[k] == 1)
      Pl_D[n++] = k;

   n_bool = Compatibility_Plasmid_Sets_Assert( Table, 
                                               Pl_R, Profile_R[Table->No_of_PLASMIDS], 
                                               Pl_D, Profile_D[Table->No_of_PLASMIDS] );
  // assert(n_bool == 1);                                              

  n_Difference_TR = 0; 
  for(k = 0; k<Table->No_of_PLASMIDS; k++)  
      if( (Profile_T[k] - Profile_R[k]) == 1 ) 
        Pl_Difference_TR[n_Difference_TR++] = k; 
  
  d = 0; 
  for(k = 0; k < n_Difference_TR; k++) 
    if (Profile_D[Pl_Difference_TR[k]] == 1) d++; 

  if( d == n_Difference_TR ) {
    /* d: Number of plasmids from the transconjuntant that are NOT in the recipient 
          but are present the donor. Effective conjugation occurs when all of them
          are transmitted (each of the them at probability Chi_C_0)

       n_0: Number of plasmids in the transconjugant that are already present in the 
            recipient. If they are present in the donor, they could enter the conjugation 
            set (with probability Chi_C_0) or not (with probability 1.0 - Chi_C_0)
        
       n: Number of plasmids in the donor that cannot go into conjugation in order to 
          produce the right transconjugant we want (each of them at probability 
          (1.0 - Chi_C). 
     */
    n_0 = Profile_T[Table->No_of_PLASMIDS] - n_Difference_TR;
    n = 0; /* No of plasmids in the Donor that are not in the Recipient
              and are none of the plasmids that should be transfered
              (those in Pl_Difference_TR[]) 
           */ 
    n_Difference_DR = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++)  
      if( (Profile_D[k] - Profile_R[k]) == 1 )
        Pl_Difference_DR[n_Difference_DR++] = k;
    
    n= 0; 
    for(k = 0; k<n_Difference_DR; k++) {
      l = 0; 
      while(l < n_Difference_TR) { 
        if(Pl_Difference_DR[k] != Pl_Difference_TR[l] )
          l++;
        else
          break; 
      }
      if(l == n_Difference_TR) n++; 
    }

    Strain_Sp_R = Table->StrainType_and_Profile[Strain_ID_R][0]; 
    Strain_Sp_T = Table->StrainType_and_Profile[Trans_ID][0];
    Strain_Sp_D = Table->StrainType_and_Profile[Strain_ID_D][0];

    assert( Strain_Sp_R == Strain_Sp_T ); 
    assert( Trans_ID > Table->n_0[Strain_Sp_T] );
    assert( Strain_ID_D > Table->n_0[Strain_Sp_D] );

    n_Invasion = 0;  
    for(k = 0; k < n_Difference_TR; k++ )
      if( Table->IBP[Strain_Sp_R][Pl_Difference_TR[k]] == 1.0 ) 
        n_Invasion++; 

    /* All plasmids in the donor the recipient requires 
       to give rise to the right transconjugant  
       should be able to invade the recipient 
       bacterial type, Strain_Sp_R.  
    */    
    assert( n_Invasion == n_Difference_TR );
    
    if ( n_Invasion != n_Difference_TR ) {
              printf(" Something is very very wrong here...\n");

              Print_Strain_Profile(Profile_T, Table->No_of_PLASMIDS); 
              printf("\t Transconjugant Profile Sp: %d ID: %d", Strain_Sp_T, Trans_ID); printf("\n");
              
              Print_Strain_Profile(Profile_R, Table->No_of_PLASMIDS);
              printf("\t Recipient Profile Sp: %d ID: %d", Strain_Sp_R, Strain_ID_R); printf("\n");

              Printf_Infection_Profile(Profile_D, Table->StrainType_and_Profile[Strain_ID_R][0], Table); 
              printf("\t\t Donor Profile (when infecting the recipient) Sp: %d ID: %d", Strain_Sp_D, Strain_ID_D); 
              printf("\n\n");

              printf("Bacteria-plasmid infection matrix: \n");
              show_DoubleMatrix(Table->IBP, Table->No_of_STRAINS, Table->No_of_PLASMIDS);
              printf("\n"); 

              for(k = 0; k < n_Difference_TR; k++ ) {
                assert( Table->IBP[Strain_Sp_R][Pl_Difference_TR[k]] == 1.0 );
                assert( Table->IBP[Strain_Sp_D][Pl_Difference_TR[k]] == 1.0 ); 
              }

              getchar();
              assert( n_Invasion == n_Difference_TR );  /* This stops the pogram and exit... */             
    }

    /* Strain_Sp_D should be able to conjugate with Strain_Sp_R */
    assert( Table->HBB[Strain_Sp_R][Strain_Sp_D] == 1.0 );  

    /* Notice that these asserts should be always passed provided that 
       both the putative recipient lists for every Strain IDs and 
       the potential donor list of every Strain ID in this putative
       recipient list have been correctly calculated. 
    */
    * Rate = pow( Table->Chi_C_0, (double)d ) * pow(1.0-Table->Chi_C_0, (double)n);    
    n_bool = 1;
  }
  else {
      * Rate = 0.0; 
      n_bool = 0;     
  }     

  free(Pl_D); free(Pl_R); free(Pl_Difference_TR);

  return (n_bool);
}                           

void Setting_Strain_Characteristic_Parameters (Parameter_Table * Table)
{
    int i, j, k, n; 
    double c, COST, RESISTANCE;  

    n = 0;
    for (i=0; i<Table->No_of_STRAINS; i++) {  
      for(j=0; j<Table->n[i]; j++) { 
        
        c = 0.0; COST = 1.0; RESISTANCE = 0.0;        
        for(k=0; k<Table->No_of_PLASMIDS; k++) {
          if(Table->Strain_Profiles[i][j][k] == 1) { 
            c = Table->Alpha_C[k];
            RESISTANCE = MAX(RESISTANCE, Table->Nu_C[k]);
            COST *= (1.0 - c);
          }
        }

        Table->Beta_AP[n]  = Table->Beta_R * COST;                                        /* Bacteria Cell Division Rates */
        Table->Eta_RP[n]   = Table->Lambda_R_1;                                           /* Bacteria Conjugation Rates   */
        Table->Delta_AP[n] = Table->Delta_R_0 + Table->Delta_R_1 * (1.0 - RESISTANCE);    /* Bacteria Death Rates         */
        Table->Mu_RP[n]    = Table->Mu;                                                   /* Bacteria Diffusion Rates     */
        Table->Segregation_Error[n]= Table->p_1;                                          /* Bacterial Segregation Error  */

        n++; 
      }
    }

    assert(n == Table->No_of_RESOURCES);
}  

void Setting_Plasmid_Characteristic_Parameters (Parameter_Table * Table)
{
    //Plasmid Cost and Resistance 
    int i; 

    for (i=0; i<Table->No_of_PLASMIDS; i++) {
      Table->Alpha_C[i]  = Table->Alpha_C_0;    /* Plasmid reproduction costs */
      Table->Nu_C[i]     = Table->Nu_C_0;       /* Plasmid resistance         */
    }
} 
