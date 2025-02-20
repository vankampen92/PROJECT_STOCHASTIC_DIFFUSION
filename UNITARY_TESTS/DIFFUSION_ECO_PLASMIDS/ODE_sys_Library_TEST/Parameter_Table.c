#include <MODEL.h>

#include "Parameter_Definitions.h"

extern gsl_rng * r; /* Global generator defined in main.c */

void  Time_Control_Upload ( Time_Control * Time, Parameter_Table * Table )
{
  /* Setup for the vector of sampling times */
  int i;

  Time->I_Time  = No_of_TIMES;

  Time->Time_0  = INITIAL_TIME;
  Time->Time_1  = FINAL_TIME;
  
  for(i=0; i < Time->I_Time; i++)
    Time->Time_Vector[i] = Time->Time_0 + (double)i * (Time->Time_1 - Time->Time_0)/(double)(Time->I_Time-1);
  
  Table->T = Time;
}

void Preparing_Initial_System_Configuration ( Parameter_Table * Table )
{
    int i, j, k, l, n, N;

    Table->TYPE_of_INITIAL_CONDITION = INITIAL_CONDITION; 
    
    Table->No_of_CELLS    = No_of_CELLS_SYSTEM; 
    
    /* Model Parameters: */ /* Model Parameters */
    Table->No_of_PROFILES = No_of_PROFILES_MAXIMUM;
    Table->No_of_PLASMIDS = No_of_PLASMIDS_MAXIMUM; 
    Table->No_of_STRAINS  = No_of_STRAINS_MAXIMUM;

    Table->K_R            = SYSTEM_SIZE; 
    Table->p_2            = SPARSITY_PARAMETER;
    Table->Delta_C_0      = COMPETITION_INDUCED_MORTALITY; 
    Table->Lambda_R_1     = COMMON_CONJUGATION_RATE;       
    Table->Beta_R         = CELL_DIVISION_RATE; 
    Table->Delta_R_0      = BASAL_DEATH_RATE; 
    Table->Delta_R_1      = STRESS_INDUCED_DEATH; 
    Table->p_1            = SEGREGATION_ERROR;
    Table->Mu             = DIFFUSION_RATE;

    Table->Chi_C_0        = PLASMID_TRASMISSION_PROBABILITY;
    Table->Alpha_C_0      = REPRODUCTION_COST;        /* Plasmid reproduction costs */
    Table->Nu_C_0         = PLASMID_RESISTANCE;       /* Plasmid resistance         */

    Table->Eta_RP    = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );
    Table->Mu_RP     = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );
    Table->Beta_AP   = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );
    Table->Delta_AP  = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );
    Table->Segregation_Error = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );
    Table->Nu_C      = (double *)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(double) );
    Table->Alpha_C   = (double *)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(double) );

    Table->Vector_Model_Variables = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double));

    /* BEGIN: Allocating and Setting up Connectivity Matrix */
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
    
      Table->Strain_Profiles = (int ***)calloc(No_of_STRAINS_MAXIMUM, sizeof(int **));
      for(j=0; j<No_of_STRAINS_MAXIMUM; j++){
        Table->Strain_Profiles[j] = (int **)calloc(No_of_PROFILES_MAXIMUM, sizeof(int *));
        for(i=0; i < No_of_PROFILES_MAXIMUM; i++) 
          Table->Strain_Profiles[j][i] = (int *)calloc(No_of_PLASMIDS_MAXIMUM + 1, sizeof(int));      
      }

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
      getchar();
      
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
      Table->MODEL_STATE_VARIABLES = Table->No_of_RESOURCES; 
      Table->LOCAL_STATE_VARIABLES = Table->No_of_RESOURCES;                                                                      

      printf(" Allocating and initializing the rest of arrays of model parameters...\n");
      Setting_Plasmid_Characteristic_Parameters (Table); //Plasmid Reproduction Costs
      Setting_Strain_Characteristic_Parameters (Table);  //Bacteria                                                                

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
      for(i=0; i<Table->No_of_PLASMIDS; i++) 
        Table->Plasmid_Compatibility_Indeces[i] = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int));     

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
    Setting_Putative_Recipient_Lists_of_Potential_Trasconjugants (Table);    
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

    printf(" Initial system configuration has been successfully established...\n");
}

void De_Allocating_Initial_System_Configuration ( Parameter_Table * Table )
{
  int i, j, k, l, n, N;
    
    free(Table->Eta_RP); 
    free(Table->Mu_RP);
    free(Table->Beta_AP);    
    free(Table->Delta_AP);
    free(Table->Segregation_Error);
    free(Table->Nu_C);
    free(Table->Alpha_C);
    
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
}      

