#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#ifndef GSL_HEADER_FILES
#define GSL_HEADER_FILES
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#endif

#define No_of_PLASMIDS_MAXIMUM 3 // 3 // 8 // 4        /* Number of Plasmids Maximum */ 

#define No_of_PROFILES_MAXIMUM 8 // 8 // 256 // 16     /* P = 2^(No_of_PLASMIDS_MAXIMUM) */

#define No_of_STRAINS_MAXIMUM 100                      /* S, number of Bacterial Strains Maximum */
                                
#define No_of_RESOURCES_MAXIMUM 800 // 800 // 25600 // 1600  /* S * P */ /* P = 2^(No_of_PLASMIDS_MAXIMUM) */
                                         /* Maximum No of Different strains and profiles (S * P), 
                                            where S is the max number of strains  
                                         */
#define SPARSITY_PARAMETER            0.5
#define COMPETITION_INDUCED_MORTALITY 0.5 
#define COMMON_CONJUGATION_RATE       1.0
#define SYSTEM_SIZE                  1000 
                                
#define VERBOSE

#define CURRENT_TIME_RDN_SEED

#include "BP_StrainsProfiles.h"

gsl_rng * r; /* Global generator defined in main.c */

int main(int argc, char **argv)
{
  int i, j, k, l, n, N;
   
  /* GNU Scientific Library:  Random numbers set up */
  const gsl_rng_type * T;
  const gsl_rng_type **t, **t0;

  gsl_rng_env_setup();  /* Environment variables  
                          GSL_RNG_TYPE and GSL_RNG_SEED are read  
                          to set the corresponding library variables 
                          gsl_rng_default and gsl_rng_default_seed. 
                        */
  T = gsl_rng_taus2;   //T = gsl_rng_default; T = gsl_rng_ranlux389;
  r = gsl_rng_alloc(T);
  t0 = gsl_rng_types_setup ();

#if defined VERBOSE
  printf("Random number information... Available generators:\n");
  for(t=t0; *t != 0; t++){
    printf("%s, ", (*t)->name);
  }
#endif
  printf("\n In this example, the random generator at work... is the '%s' generator\n\n",
         gsl_rng_name(r));
#if defined CURRENT_TIME_RDN_SEED
  GSL_Init_Random_Seed(r);
#else
  GSL_Init_Random_Seed_from_File(r);
#endif

  /* End of Random number setting */

  #if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); //getchar();
  /*   END: Checking Random Number Generator Setup */
  #endif

  /* Allocating and Setting Matrices:
            1.  Competition Matrix, ABB 
            2.  Conjugation Matrix, HBB 
            3.  Plasmid-Bacteria Infection Matrix, IPB
            4.  Plasmid-Plasmid Compatibility Matrix, CPP
  */
    Parameter_Table * Table = (Parameter_Table *)calloc(1, sizeof(Parameter_Table));
    /* n_0[i]: Strain ID of the plasmid-free profile 
               associated to the i-th bacterial type */
    Table->n_0 = (int *)calloc(No_of_STRAINS_MAXIMUM, sizeof(int));
    
    /* n_R[i]: No of Recipients of Strain ID 'i' */ 
    Table->n_R = (int *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(int));
    
    /* No of actual viable profiles per strain. It can be different per each strain */ 
    Table->n = (int *)calloc(No_of_STRAINS_MAXIMUM, sizeof(int));

    // Table->StrainType_and_Profile[N][0] = i;     /* Strain Type          */
    // Table->StrainType_and_Profile[N][1] = n[i];  /* Profile No Specifier */
    Table->StrainType_and_Profile = (int **)calloc(No_of_RESOURCES_MAXIMUM, sizeof(int *));
    for(j=0; j<No_of_RESOURCES_MAXIMUM; j++)
      Table->StrainType_and_Profile[j] = (int *)calloc(2, sizeof(int));     

    Table->Competition_Induced_Death = (double **)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double *));
    for(j=0; j<No_of_RESOURCES_MAXIMUM; j++)
      Table->Competition_Induced_Death[j] = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double));    

    Table->ABB = (double **)calloc(No_of_STRAINS_MAXIMUM, sizeof(double *) ); 
    for(i=0; i<No_of_STRAINS_MAXIMUM; i++) 
      Table->ABB[i] = (double *)calloc(No_of_STRAINS_MAXIMUM, sizeof(double) );
    
    Table->HBB = (double **)calloc(No_of_STRAINS_MAXIMUM, sizeof(double *) ); 
    for(i=0; i<No_of_STRAINS_MAXIMUM; i++) 
      Table->HBB[i] = (double *)calloc(No_of_STRAINS_MAXIMUM, sizeof(double) );

    Table->IBP = (double **)calloc(No_of_STRAINS_MAXIMUM, sizeof(double *) ); 
    for(i=0; i<No_of_STRAINS_MAXIMUM; i++) 
      Table->IBP[i] = (double *)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(double) );
    
    Table->CPP = (double **)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(double *) ); 
    for(i=0; i<No_of_PLASMIDS_MAXIMUM; i++) 
      Table->CPP[i] = (double *)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(double) );
   /* END: Allocating and Setting Matrices */ 

    Table->Strain_Profiles = (int ***)calloc(No_of_STRAINS_MAXIMUM, sizeof(int **));
    for(j=0; j<No_of_STRAINS_MAXIMUM; j++){
      Table->Strain_Profiles[j] = (int **)calloc(No_of_PROFILES_MAXIMUM, sizeof(int *));
      for(i=0; i < No_of_PROFILES_MAXIMUM; i++) 
        Table->Strain_Profiles[j][i] = (int *)calloc(No_of_PLASMIDS_MAXIMUM + 1, sizeof(int));      
    }

    Table->No_of_PROFILES = No_of_PROFILES_MAXIMUM;
    Table->No_of_PLASMIDS = No_of_PLASMIDS_MAXIMUM; 
    Table->No_of_STRAINS  = No_of_STRAINS_MAXIMUM;
    Table->p_2            = SPARSITY_PARAMETER;
    Table->Delta_C_0      = COMPETITION_INDUCED_MORTALITY; 
    Table->K_R            = SYSTEM_SIZE; 
    Table->Lambda_R_1     = COMMON_CONJUGATION_RATE; 

    Table->Eta_RP = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));

    printf(" Parameter_Table structure has been correctly allocated and initiated\n");

    Setting_Interaction_Matrices (Table);

    /* B E G I N :  Writing out interaction matrices: */
    printf("Competition matrix: \n");
    show_DoubleMatrix(Table->ABB, Table->No_of_STRAINS, Table->No_of_STRAINS);
    printf("\n");

    printf("Conjugation matrix: \n");
    show_DoubleMatrix(Table->HBB, Table->No_of_STRAINS, Table->No_of_STRAINS);
    printf("\n");

    printf("Bacteria-plasmid infection matrix: \n");
    show_DoubleMatrix(Table->IBP, Table->No_of_STRAINS, Table->No_of_PLASMIDS);
    printf("\n");    

    printf("Plasmid-plasmid compatibilty matrix: \n");
    show_DoubleMatrix(Table->CPP, Table->No_of_PLASMIDS, Table->No_of_PLASMIDS);
    printf("\n");

    /* Determining actual No_of_RESOURCES after considering the constraints established 
       by the plasmid-plasmid compatibilty matrix (certain plasmids are incompatible in the same 
       bacterial cell), and infection matrix (certain strains are immune to certain plasmids) 
    */
    Table->No_of_RESOURCES = Determining_actual_No_of_RESOURCES (Table); /* No_of_LOCAL_VARIABLES, 
                                                                            this is, for instance, 
                                                                            the i-th strain with 
                                                                            profile p 
                                                                          */
                                                                          /* In Ying-Jie notation, 
                                                                             No of SUBPOPULATIONS 
                                                                          */
    Table->Eta_RP = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));
    for(i=0; i<Table->No_of_RESOURCES; i++) 
      Table->Eta_RP[i] = Table->Lambda_R_1;

    /* Both allocating and initializing Adjacancy Lists from interaction matrices */
    Table->Competition_List_Indeces = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
    Table->Conjugation_List_Indeces = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
    Table->Recipient_List_Indeces   = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
    Table->Putative_Recipient_List_Indeces = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
    Table->Donor_List_Indeces       = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
    for(i=0; i<Table->No_of_RESOURCES; i++) {
      Table->Competition_List_Indeces[i] = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
      Table->Conjugation_List_Indeces[i] = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
      Table->Recipient_List_Indeces[i]   = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
      Table->Putative_Recipient_List_Indeces[i]   = (int *)calloc(Table->No_of_PROFILES + 1, sizeof(int));
      Table->Donor_List_Indeces[i]       = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
    }

    Table->Plasmid_Compatibility_Indeces = (int **)calloc(Table->No_of_PLASMIDS, sizeof(int *));
    for(i=0; i<Table->No_of_PLASMIDS; i++) {
      Table->Plasmid_Compatibility_Indeces[i] = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int)); 
    }

    printf(" Initialization of the adjancy lists from the interaction matrices... (pres key)"); getchar();
    Setting_Adjacency_Lists_from_Interaction_Matrices (Table);
    printf(" Success!!! List initializations concluded!!!\n");
    printf("\n");

    printf("Printing Strains, Profiles and Lists:\n"); 
    // getchar();

    Printing_Strains_Profiles_and_Lists(Table);
  
    printf("\n");
    printf("Calculating Putatitive Recipient Lists of Potential Transconjugants...\n");
    getchar();
    Setting_Putatitive_Recipient_Lists_of_Potential_Trasconjugants (Table);    
    getchar();
    Printing_Putative_Recipient_Lists(Table);

    Table->DoRe = (Donor_Recipient_Pair ****)calloc(Table->No_of_RESOURCES, sizeof(Donor_Recipient_Pair ***));
    for(k = 0; k < Table->No_of_RESOURCES; k++)
      Table->DoRe[k] = (Donor_Recipient_Pair ***)calloc(Table->Putative_Recipient_List_Indeces[k][Table->No_of_PROFILES], 
                                                        sizeof(Donor_Recipient_Pair **));    
    /* The rest of required allocation and initialization is done within the function:

       Setting_Reactive_Recipient_Donor_Pairs_and_Rates (...);       
    */  
    printf(" Calculating the set of effective putative recipient-donor pairs for every transconjugant Strain ID...\n"); 
    printf(" (notice that plasmid-free Strain IDs can never result from conjugation)\n"); getchar(); 
    Setting_Reactive_Recipient_Donor_Pairs_and_Rates(Table);
      
    getchar();
    Printing_Putative_Recipient_Lists(Table);

    getchar(); 
    printf("Bacteria-plasmid infection matrix: \n");
    show_DoubleMatrix(Table->IBP, Table->No_of_STRAINS, Table->No_of_PLASMIDS);
    printf("\n");    

    printf("Plasmid-plasmid compatibilty matrix: \n");
    show_DoubleMatrix(Table->CPP, Table->No_of_PLASMIDS, Table->No_of_PLASMIDS);
    printf("\n");

    /* De-Allocating Matrices and lists:
            1.  Competition Matrix, ABB 
            2.  Conjugation Matrix, HBB 
            3.  Plasmid-Bacteria Infection Matrix, IPB
            4.  Plasmid-Plasmid Compatibility Matrix, CPP
    */
    for(i=0; i<No_of_STRAINS_MAXIMUM; i++) {
      free(Table->ABB[i]); 
      free(Table->HBB[i]);
      free(Table->IBP[i]); 
    }
    free(Table->ABB);  
    free(Table->HBB);  
    free(Table->IBP); 

    for(i=0; i<No_of_PLASMIDS_MAXIMUM; i++) 
      free(Table->CPP[i]);
    free(Table->CPP);  

    /* De-Allocating the whole Table->DoRe 'Donor_Recipient_Pair' data structure */
    for(k = 0; k < Table->No_of_RESOURCES; k++) {
      for(l = 0; l < Table->Putative_Recipient_List_Indeces[k][Table->No_of_PROFILES] ; l++) {
        
        int Strain_ID_R = Table->Putative_Recipient_List_Indeces[k][l];
        printf(" Number of effective donors of the %d-th putative recipient of strain ID %d:\t %d (out of %d)\n", l, k,  
                 Table->DoRe[k][l][0]->N, 
                 Table->Donor_List_Indeces[Strain_ID_R][Table->No_of_RESOURCES]); 

        n == 0;    N = Table->DoRe[k][l][0]->N;  
        while( n < N ) free(Table->DoRe[k][l][n++]);

        free(Table->DoRe[k][l]);   
      }
      free(Table->DoRe[k]);
    }     
    free(Table->DoRe);

    /* De-Allocating Strain Profiles and lists */
    for(j=0; j<No_of_STRAINS_MAXIMUM; j++){

      for(i=0; i < No_of_PROFILES_MAXIMUM; i++) 
        free(Table->Strain_Profiles[j][i]);      
      free(Table->Strain_Profiles[j]);
    }
    free(Table->Strain_Profiles);
    
    for(i=0; i<Table->No_of_RESOURCES; i++) {
      free(Table->Competition_List_Indeces[i]);
      free(Table->Conjugation_List_Indeces[i]);
      free(Table->Recipient_List_Indeces[i])  ;
      free(Table->Putative_Recipient_List_Indeces[i])  ;
      free(Table->Donor_List_Indeces[i])      ;
    }    
    free(Table->Competition_List_Indeces) ;
    free(Table->Conjugation_List_Indeces) ;
    free(Table->Recipient_List_Indeces)   ;
    free(Table->Putative_Recipient_List_Indeces)   ;
    free(Table->Donor_List_Indeces)       ;

    free(Table->Eta_RP);

    for(i=0; i<Table->No_of_PLASMIDS; i++) 
      free(Table->Plasmid_Compatibility_Indeces[i]); 
    free(Table->Plasmid_Compatibility_Indeces);

    for(j=0; j<No_of_RESOURCES_MAXIMUM; j++)
      free(Table->StrainType_and_Profile[j]);       
    free(Table->StrainType_and_Profile);

    for(j=0; j<No_of_RESOURCES_MAXIMUM; j++)
      free(Table->Competition_Induced_Death[j]);
    free(Table->Competition_Induced_Death);

    free(Table->n);    /* No of actual viable profiles per 
                          bacterial type, different for each of them. */
    free(Table->n_0);  /* Strain ID of the plasmid-fre profiles 
                          for each bacterial type */
    free(Table->n_R);  /* No of potential recipients per each bacterial 
                          type. This is redundant to the same number 
                          stored in:

                          Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES], 
                          
                          which also stores the length of the recipient list 
                          for the i-th Strain ID. 
                        */

    printf("\n ------------------------------------------------------------------------------ \n");
    printf("No of RESOURCES (No of STRAIN IDs): %d\n", Table->No_of_RESOURCES);
    printf("No of STRAIN TYPES: %d\n", Table->No_of_STRAINS);
    printf("No of PLASMIDS: %d\n", Table->No_of_PLASMIDS);
    printf("No of TOTAL EVENTS: %d\n", Table->No_of_STRAINS * 6 + Table->No_of_CONJUGATION_EVENTS);
    printf("No of Conjugation Pairs (Donor, Recipient): %d\n", Table->No_of_CONJUGATION_EVENTS);

    free(Table); 
    
    printf(" ---> End of Program\n");
    return(0);
}


