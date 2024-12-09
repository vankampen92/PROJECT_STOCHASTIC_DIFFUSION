 #include <MODEL.h>

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