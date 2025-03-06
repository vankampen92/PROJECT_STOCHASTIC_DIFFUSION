#include <MODEL.h>

#if defined CPGPLOT_REPRESENTATION
#include <include.CPG.extern.h>
#endif

#include <include.Parameter_Model.extern.h>
#include <include.Output_Variables.extern.h>
#include <include.Parameter_Space.extern.h>
#include <include.Time_Control.extern.h>
#include <include.Initial_Conditions.extern.h>
#include <include.Error_Control.extern.h>
#include <include.Type_of_Model.extern.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C( Parameter_Table * Table )
{
  int i, j, a;
  int MODEL_OUTPUT_VARIABLES;   /* Actual no of MODEL (output) VARIABLES */
  int MODEL_INPUT_PARAMETERS;  /* Actual no of MODEL (input) PARAMETERS */
  int MODEL_PARAMETERS_EXTENDED; 

  MODEL_INPUT_PARAMETERS = MODEL_PARAMETERS_MAXIMUM;
  MODEL_OUTPUT_VARIABLES = OUTPUT_VARIABLES_MAXIMUM;
  MODEL_PARAMETERS_EXTENDED = MODEL_PARAMETERS_MAXIMUM * (1 + No_of_TDC_FUNC_AUX_PARAM_MAX); 
  // int MODEL_STATE_VARIABLES;  /* Actual no of MODEL (state) VARIABLES 
  /* Model parameters are all input parameters defining and controling
     model dynamics (both through simulation and mathematically, in case
     the model has a mathematical description such as a system of ODEs)
  */
  Table->Code_Parameters = (char **)calloc( MODEL_PARAMETERS_EXTENDED, sizeof(char *) );
  Table->Symbol_Parameters = (char **)calloc(MODEL_PARAMETERS_EXTENDED, sizeof(char *) );
  Table->Symbol_Parameters_Greek_CPGPLOT = (char **)calloc( MODEL_PARAMETERS_EXTENDED,
							    sizeof(char *));
  Table->Name_Parameters = (char **)calloc( MODEL_PARAMETERS_EXTENDED, sizeof(char *) );
  for (i=0; i<MODEL_PARAMETERS_EXTENDED; i++){
    Table->Code_Parameters[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Name_Parameters[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Symbol_Parameters[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Symbol_Parameters_Greek_CPGPLOT[i] = (char *)malloc( 100 * sizeof(char) );
  }
  Table->Default_Vector_Parameters = (double *)calloc(MODEL_PARAMETERS_EXTENDED, sizeof(double) );
  Table->Vector_Parameters = (double *)calloc( MODEL_PARAMETERS_EXTENDED, sizeof(double) );
  Table->Index = (int *)calloc( MODEL_PARAMETERS_EXTENDED, sizeof(int) );

  /* Output Variables are any measure of the state model variables of any function of
     these at any given time:
     MODEL_OUTPUT_VARIABLES is the Total Number of Total Output Variables
     (as it appears in defintion_OutPut_Variables.c file)( )
   */
  Table->Output_Variable_Name = (char **)malloc( MODEL_OUTPUT_VARIABLES * sizeof(char *) );
  Table->Output_Variable_Symbol = (char **)malloc( MODEL_OUTPUT_VARIABLES * sizeof(char *) );
  for (i=0; i<MODEL_OUTPUT_VARIABLES; i++){
    Table->Output_Variable_Name[i] = (char *)malloc( 100 * sizeof(char) );
    Table->Output_Variable_Symbol[i] = (char *)malloc( 100 * sizeof(char) );
  }
  Table->OUTPUT_VARIABLE_INDEX = (int *)malloc( SUB_OUTPUT_VARIABLES * sizeof(int) );
  Table->Vector_Output_Variables = (double *)malloc( SUB_OUTPUT_VARIABLES * sizeof(double) );
  Table->Matrix_Output_Variables = (double **)malloc( SUB_OUTPUT_VARIABLES * sizeof(double *) );
  for (i=0; i<SUB_OUTPUT_VARIABLES; i++)
    Table->Matrix_Output_Variables[i] = (double *)malloc( I_Time * sizeof(double) );

  Table->Default_Vector_Output_Variables = (double *)malloc( MODEL_OUTPUT_VARIABLES * sizeof(double) );

  Table->Vector_Model_Variables_MultiStability =(double **)calloc( 3, sizeof(double *) );

  /* Model Variables represent only the state model variables, i.e., the set of variables
     completely defining the state of the system at any given time. Model Variables can be,
     of course, also output variables. They can be defined as output variables
     whose definition can be found at Definition_Output_Variables.c */

  /* Some implementations of this code require to alloc memmory according to
     a number of state variables that can change dynamically. In that case, this
     6 lines of code should be commented out. Notice alse the corresponding lines
     in void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * T )

     For most calculations, I opted for allocating this part of the data structure
     only if needed (see ./MODEL_CALCULATIONS/TIME_EVO_DETERMINISTIC/MODEL.c) regardless
     the number of MODEL STATE VARIABLES changing or not during execution. Essentially,
     the numerical integration of a system of ODEs relies on Table->Vector_Model_Variables.
     Therefore, a driver, usaully called void MODEL(Parameter_Table *), is in charge
     of doing this allocation previous to calling the numerical integration function,
     and deallocating this memmory right after.
  */
  // Table->Model_Variable_Name = (char **)malloc( MODEL_STATE_VARIABLES * sizeof(char *) );
  // for (i=0; i<MODEL_STATE_VARIABLES; i++){
  //   Table->Model_Variable_Name[i] = (char *)malloc( 100 * sizeof(char) );
  // }
  // Table->Vector_Model_Variables = (double *)malloc( MODEL_STATE_VARIABLES * sizeof(double) );
  // Table->Vector_Model_Variables_Stationarity = (double *)malloc( MODEL_STATE_VARIABLES * sizeof(double) );

  if(     TYPE_of_NETWORK == 0) No_of_NEIGHBORS = No_of_CELLS - 1;
  else if(TYPE_of_NETWORK == 1) No_of_NEIGHBORS = 4;            /* Von Neumann */
  else {
    printf(" No more network types have been defined so far\n");
    printf(" The program will exit\n");
    exit(0); 
  }
  
  Table->Lambda_R  = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) ); /* DIFFUSION_S_RESOURCES   */
  Table->Delta_R   = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) ); /* DIFFUSION_S_RESOURCES   */
  Table->Delta_RP  = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) ); /* DIFFUSION_ECOEVO_PLANTS */

  Table->Alpha_C   = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );
  Table->Nu_C   = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );

  #if defined(DIFFUSION_ECO_PLASMIDS) || defined(DIFFUSION_ECO_1B1P)
  assert( No_of_RESOURCES_MAXIMUM > No_of_PLASMIDS_MAXIMUM );
  #endif 

  Table->Nu_Consumers    = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) ); /* DIFFUSION_HII_nD */
  Table->Theta_Consumers = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) ); /* DIFFUSION_HII_nD */
  Table->y_R_i     = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );       /* DIFFUSION_HII_nD */
  
  Table->Segregation_Error = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) ); /* DIFFUSION_ECO_PLASMIDS and _ECO_1B1P */
  Table->Eta_RP    = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );         /* DIFFUSION_ECO_PLASMIDS and _ECO_1B1P */   
  Table->Mu_RP     = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );         /* DIFFUSION_ECO_PLASMIDS and _ECO_1B1P */
  Table->Beta_AP   = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );         /* DIFFUSION_ECO_PLASMIDS and _ECO_1B1P */ 
  Table->Delta_AP  = (double *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(double) );         /* DIFFUSION_ECO_PLASMIDS and _ECO_1B1P */

  /* BEGIN: Allocating and Setting up Connectivity Matrix */
  Table->Metapop_Connectivity_Matrix = (double ***)calloc(No_of_RESOURCES_MAXIMUM,
							  sizeof(double **) );
  for(a=0; a<No_of_RESOURCES_MAXIMUM; a++) { 
    Table->Metapop_Connectivity_Matrix[a] = (double **)calloc(No_of_CELLS,
							      sizeof(double *) );
    for(i=0; i<No_of_CELLS; i++) {
      Table->Metapop_Connectivity_Matrix[a][i] = (double *)calloc(No_of_NEIGHBORS+1,
								  sizeof(double) );
    }
  }
  /* END: Allocating and Setting up Connectivity Matrix */

#if defined(DIFFUSION_ECO_PLASMIDS) || defined(DIFFUSION_ECO_1B1P)
  #if defined DIFFUSION_ECO_PLASMIDS
  /* BEGIN: Allocating and Setting Matrices:
            1.  Competition Matrix, ABB 
            2.  Conjugation Matrix, HBB 
            3.  Plasmid-Bacteria Infection Matrix, IPB
            4.  Plasmid-Plasmid Compatibility Matrix, CPP
  */
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
  #endif
    /* No of actual viable profiles per strain. It can be different per each strain */ 
    Table->n = (int *)calloc(No_of_STRAINS_MAXIMUM, sizeof(int)); 
    Table->n_0 = (int *)calloc(No_of_STRAINS_MAXIMUM, sizeof(int)); /* Strain ID corresponding 
                                                                       to a strain without plasmids 
                                                                    */
    Table->n_R = (int *)calloc(No_of_RESOURCES_MAXIMUM, sizeof(int)); /* No of Recipeints of 
                                                                         every Strain ID 
                                                                      */
    Table->Strain_Profiles = (int ***)calloc(No_of_STRAINS_MAXIMUM, sizeof(int **));
    for(j=0; j<No_of_STRAINS_MAXIMUM; j++){
      Table->Strain_Profiles[j] = (int **)calloc(No_of_PROFILES_MAXIMUM, sizeof(int *));
      for(i=0; i < No_of_PROFILES_MAXIMUM; i++) 
        Table->Strain_Profiles[j][i] = (int *)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(int));      
    }

    Table->Global_Strains_Population = (double *)calloc(No_of_STRAINS_MAXIMUM, sizeof(double));  
    Table->Global_Plasmid_Population = (double *)calloc(No_of_PLASMIDS_MAXIMUM, sizeof(double));  
#endif  
  
  Table->A_P = (int *)calloc(N_E, sizeof(int) );
  Table->RA_P = (int *)calloc(N_E, sizeof(int) );
  Table->ARA_P = (int **)calloc(N_E, sizeof(int *) );
  for(i=0; i<N_E; i++)
    Table->ARA_P[i] = (int *)calloc(N_E, sizeof(int) );
}

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * Table )
{
  int i, j, a;

  for(i=0; i<Table->N_E; i++)
    free(Table->ARA_P[i]);

  free(Table->A_P);
  free(Table->RA_P);
  free(Table->ARA_P);

  for (i=0; i < Table->MODEL_INPUT_PARAMETERS; i++){
    free( Table->Code_Parameters[i] );
    free( Table->Name_Parameters[i] );
    free( Table->Symbol_Parameters[i] ); 
    free( Table->Symbol_Parameters_Greek_CPGPLOT[i] ); 
  }
  free(Table->Code_Parameters);
  free(Table->Name_Parameters);
  free(Table->Symbol_Parameters); 
  free(Table->Symbol_Parameters_Greek_CPGPLOT); 
  
  free(Table->Default_Vector_Parameters);

  free(Table->Index);
  free(Table->Vector_Parameters);

  for (i=0; i < OUTPUT_VARIABLES_MAXIMUM; i++){
    free(Table->Output_Variable_Name[i]);
    free(Table->Output_Variable_Symbol[i]);
  }
  free(Table->Output_Variable_Symbol);
  free(Table->Output_Variable_Name);

  for (i=0; i<Table->MODEL_STATE_VARIABLES; i++){
     free(Table->Model_Variable_Name[i]);
     free(Table->Model_Variable_Symbol[i]);
  }
  free(Table->Model_Variable_Name);
  free(Table->Model_Variable_Symbol);

  /* Some implementations of this code require to alloc memmory according to
     a number of model variables that can change dynamically.

     For most calculations,  I opted for allocating/de-allocating this part of
     the data structure only if needed
     (see, for instance, ./MODEL_CALCULATIONS/TIME_EVO_DETERMINISTIC/MODEL.c)
     regardless the number of MODEL STATE VARIABLES changing or not during execution.
  */

  // for (i=0; i<MODEL_STATE_VARIABLES; i++){
  //   free(Table->Model_Variable_Name[i]);
  // }
  // free(Table->Model_Variable_Name);
  // free( Table->Vector_Model_Variables );
  // free( Table->Vector_Model_Variables_Stationarity );
  // free( Table->Vector_Model_Variables_MultiStability[0] );
  // free( Table->Vector_Model_Variables_MultiStability[1] );
  // free( Table->Vector_Model_Variables_MultiStability[2] );
  free(Table->Vector_Model_Variables_MultiStability);

  free(Table->OUTPUT_VARIABLE_INDEX);
  free(Table->Vector_Output_Variables);
  for (i=0; i<Table->SUB_OUTPUT_VARIABLES; i++) 
    free(Table->Matrix_Output_Variables[i]);
  free(Table->Matrix_Output_Variables);

  free(Table->Default_Vector_Output_Variables);

  free(Table->Lambda_R);
  free(Table->Delta_R);

  free(Table->Alpha_C);
  free(Table->Nu_C);

  free(Table->Nu_Consumers);
  free(Table->Theta_Consumers);
  free(Table->y_R_i);

  free(Table->Segregation_Error);
  free(Table->Delta_RP); 
  free(Table->Eta_RP);
  free(Table->Mu_RP);
  free(Table->Beta_AP);
  free(Table->Delta_AP);

  for(a=0; a<No_of_RESOURCES_MAXIMUM; a++) {  
    for(i=0; i<Table->No_of_CELLS; i++) 
      free(Table->Metapop_Connectivity_Matrix[a][i]); 
    
    free(Table->Metapop_Connectivity_Matrix[a]); 
  }
  free(Table->Metapop_Connectivity_Matrix);

  #if defined(DIFFUSION_ECO_PLASMIDS) || defined(DIFFUSION_ECO_1B1P)
    #if defined DIFFUSION_ECO_PLASMIDS
    /* BEGIN: De-Allocating Matrices:
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
    /* END: Allocating and Setting up Connectivity Matrix */
    #endif

  /* No of actual viable profiles per strain. It can be different per each strain */ 
    free(Table->n);  
    free(Table->n_0);
    free(Table->n_R);
   
    for(j=0; j<No_of_STRAINS_MAXIMUM; j++){
      for(i=0; i < No_of_PROFILES_MAXIMUM; i++) 
        free(Table->Strain_Profiles[j][i]);      

      free(Table->Strain_Profiles[j]);
    }
    free(Table->Strain_Profiles);

    free(Table->Global_Strains_Population);  
    free(Table->Global_Plasmid_Population);
  #endif
  
  #if defined DIFFUSION_ECO_PLASMIDS      
    for(i=0; i<Table->No_of_RESOURCES; i++) {
      free(Table->Competition_Induced_Death[i]);
      free(Table->Competition_List_Indeces[i]);
      free(Table->Conjugation_List_Indeces[i]);
      free(Table->Recipient_List_Indeces[i])  ;
      free(Table->Donor_List_Indeces[i])      ;
      free(Table->No_of_Event_Conjugation_Pair[i]);
      free(Table->StrainType_and_Profile[i]);
    }
    free(Table->StrainType_and_Profile);
    free(Table->No_of_Event_Conjugation_Pair);
    free(Table->Competition_Induced_Death);
    free(Table->Competition_List_Indeces) ;
    free(Table->Conjugation_List_Indeces) ;
    free(Table->Recipient_List_Indeces)   ;
    free(Table->Donor_List_Indeces)       ;

    for(i=0; i<Table->No_of_PLASMIDS; i++) 
      free(Table->Plasmid_Compatibility_Indeces[i]); 
    free(Table->Plasmid_Compatibility_Indeces);

    for(i=0; i<Table->No_of_CONJUGATION_EVENTS; i++) 
      free(Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[i]); 
    free(Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs);    
  #endif  
}

void P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( Parameter_Table * Table, int * Index_Output_Variables )
{
  int i, j, a;
  
  /* Stochastic Realizations */
  Table->Realizations      = Realizations;

  Table->No_of_CELLS       = No_of_CELLS; 
  Table->No_of_INDIVIDUALS = No_of_INDIVIDUALS;
  Table->No_of_CELLS_X     = No_of_CELLS_X;
  Table->No_of_CELLS_Y     = No_of_CELLS_Y;

  /* Type of Model upload  */
  Table->TYPE_of_MODEL = TYPE_of_MODEL;  
  assert_right_model_definition( Table );

  /* Parameter Model Upload */
  Parameter_Values_into_Parameter_Table(Table);

  #if defined(DIFFUSION_ECO_PLASMIDS) || defined(DIFFUSION_ECO_1B1P)
    assert(Table->TYPE_of_MODEL == 21 || Table->TYPE_of_MODEL == 22); 

    #if defined DIFFUSION_ECO_PLASMIDS  
    /* Setting up of the strain and plasmid characteristics and the 4 driving interaction matrices 
     The sparcity of these matrix determine the number of state variables to describe the configuration 
     of the system */
      assert(Table->TYPE_of_MODEL == 21); 

      Setting_Interaction_Matrices (Table); // ABB, HBB, IBP, CPP
  
    /* Determining actual No_of_RESOURCES after considering the constraints established 
       by the plasmid-plasmid compatibilty matrix (certain plasmids are incompatible in the same 
       bacterial cell), and infection matrix (certain strains are immune to certain plasmids) 
    */
      Table->No_of_RESOURCES = Determining_actual_No_of_RESOURCES (Table); /* No_of_LOCAL_VARIABLES, 
                                                                            this is, for instance, 
                                                                            the i-th strain with 
                                                                            profile p 
                                                                          */
                                                                          /* In Ying-Jie labelling, 
                                                                             No of SUBPOPULATIONS 
                                                                          */
    #endif
    #if defined DIFFUSION_ECO_1B1P
      assert(Table->TYPE_of_MODEL == 22); 

      Table->n[0]   = 2;  
      Table->n_0[0] = 0;     
      Table->n_R[0] = 0; Table->n_R[1] = 1;  /* No of Recipeints of every Strain ID */

      for(j=0; j<Table->No_of_STRAINS; j++){
          Table->Strain_Profiles[j][0][0] = 0; 
          Table->Strain_Profiles[j][1][0] = 1;       
      }
      
      Table->No_of_CONJUGATION_EVENTS = 1;
    #endif

    Setting_Plasmid_Characteristic_Parameters (Table); //Plasmid Reproduction Costs
    Setting_Strain_Characteristic_Parameters (Table);  //Bacteria

    #if defined DIFFUSION_ECO_PLASMIDS
      assert(Table->TYPE_of_MODEL == 21);

      /* Both allocating and initializing Adjacancy Lists from interaction matrices */
      Table->Competition_Induced_Death = (double **)calloc(Table->No_of_RESOURCES, sizeof(double *));
      Table->Competition_List_Indeces = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
      Table->Conjugation_List_Indeces = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
      Table->Recipient_List_Indeces   = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
      Table->Donor_List_Indeces       = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
      Table->Plasmid_Compatibility_Indeces = (int **)calloc(Table->No_of_PLASMIDS, sizeof(int *));
      Table->No_of_Event_Conjugation_Pair = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
      Table->StrainType_and_Profile       = (int **)calloc(Table->No_of_RESOURCES, sizeof(int *));
      for(i=0; i<Table->No_of_RESOURCES; i++) {
        Table->StrainType_and_Profile[i]       = (int *)calloc(2, sizeof(int));
        Table->No_of_Event_Conjugation_Pair[i] = (int *)calloc(Table->No_of_RESOURCES, sizeof(int));
        Table->Competition_Induced_Death[i]    = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));
        Table->Competition_List_Indeces[i]     = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
        Table->Conjugation_List_Indeces[i]     = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
        Table->Recipient_List_Indeces[i]       = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
        Table->Donor_List_Indeces[i]           = (int *)calloc(Table->No_of_RESOURCES + 1, sizeof(int));
      }
      for(i=0; i<Table->No_of_PLASMIDS; i++) 
        Table->Plasmid_Compatibility_Indeces[i] = (int *)calloc(Table->No_of_PLASMIDS + 1, sizeof(int)); 
        
      Setting_Adjacency_Lists_from_Interaction_Matrices (Table); /* This sets up Table->No_of_CONJUGATION_EVENTS */

      Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs = (int **)calloc(Table->No_of_CONJUGATION_EVENTS, 
                                                                              sizeof(int *));
      for(i=0; i<Table->No_of_CONJUGATION_EVENTS; i++) 
        Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[i] = (int *)calloc(2, sizeof(int));
    #endif 
  #endif

  Model_Variables_Code_into_Parameter_Table (Table);

  /* Total number of potential input paramters */
  Table->MODEL_INPUT_PARAMETERS   = MODEL_PARAMETERS_MAXIMUM;
  Table->OUTPUT_VARIABLES_GENUINE = Table->LOCAL_STATE_VARIABLES + OUTPUT_VARIABLES_TRUE_DERIVED;
  /* Three are the different cases in definition_OutPut_Variables.c */

  /* Total number of potential state variables */
  /* Total number of potential state variables */
  Table->MODEL_STATE_VARIABLES = Table->K + 1;
  /* Total number of potential output variables */
  Table->MODEL_OUTPUT_VARIABLES = Table->OUTPUT_VARIABLES_GENUINE + Table->MODEL_STATE_VARIABLES;

  /* Total number of actual model output variables */
  Table->SUB_OUTPUT_VARIABLES   = SUB_OUTPUT_VARIABLES;
  
  /// Setting MODEL INPUT PARAMETERS up !!!
  Table->Growth_Function_Type = Growth_Function_Type;
  /* BEGIN: Parameter default values into vector structure */
  for(i = 0; i<Table->MODEL_INPUT_PARAMETERS; i++){
    Table->Default_Vector_Parameters[i] = AssignStructValue_to_VectorEntry(i, Table);
  }
  /*   END: Parameter default values into vector structure */

  /* BEGIN: Names and codes for model parameters */
  for(i = 0; i<Table->MODEL_INPUT_PARAMETERS; i++){
    AssignLabel_to_Model_Parameters(i, Table->Name_Parameters[i], Table);
    AssignCodes_to_Model_Parameters(i, Table->Code_Parameters[i], Table);
    AssignSymbol_to_Model_Parameters(i, Table->Symbol_Parameters[i], Table);
    AssignCPGPLOT_Symbol_to_Model_Parameters(i, Table->Symbol_Parameters_Greek_CPGPLOT[i], Table); 
  }
  /*   END: Names and codes for model parameter  */

  /* BEGIN: Names for aux model parameters (time dependence) */
  int i_0 = Table->MODEL_INPUT_PARAMETERS; 
  int i_1 = Table->MODEL_INPUT_PARAMETERS * (1 + No_of_TDC_FUNC_AUX_PARAM_MAX); 
  for(i = i_0; i < i_1; i++){
    AssignLabel_to_Model_Parameters(i, Table->Name_Parameters[i], Table);
    /* AssignCodes_to_Model_Parameters(i, Table->Code_Parameters[i], Table); */
    AssignSymbol_to_Model_Parameters(i, Table->Symbol_Parameters[i], Table);
    AssignCPGPLOT_Symbol_to_Model_Parameters(i, Table->Symbol_Parameters_Greek_CPGPLOT[i],Table); 
  }
  /*   END: Names and codes for model parameter  */
  
  /// Setting MODEL STATE VARIABLES up !!!  
  int MODEL_STATE_VARIABLES    = Table->K + 1;
  Table->Model_Variable_Name   = (char **)malloc( MODEL_STATE_VARIABLES * sizeof(char *) );
  Table->Model_Variable_Symbol = (char **)malloc( MODEL_STATE_VARIABLES * sizeof(char *) );
  for (i=0; i<MODEL_STATE_VARIABLES; i++){
     Table->Model_Variable_Name[i] = (char *)malloc( 100 * sizeof(char) );
     Table->Model_Variable_Symbol[i] = (char *)malloc( 100 * sizeof(char) );
  }
  /* BEGIN: Names and symbols for model state variables */
  for(i = 0; i<Table->MODEL_STATE_VARIABLES; i++){
    AssignLabel_to_Model_Variables(i, Table->Model_Variable_Name[i], Table);
    AssignSymbol_to_Model_Variables(i, Table->Model_Variable_Symbol[i], Table);
  }
  /*   END: Names and codes for model variables  */

  /// Setting MODEL OUTPUT VARIABLES up !!!  
  /* BEGIN: Names and symbol for output variables           */
  for(i = 0; i<Table->MODEL_OUTPUT_VARIABLES; i++){
    AssignLabel_to_Output_Variables(i, Table->Output_Variable_Name[i], Table);
    // AssignSymbol_to_Output_Variables(i, Table->Output_Variable_Symbol[i], Table);
    AssignCPGPLOT_Symbol_to_Output_Variables(i, Table->Output_Variable_Symbol[i], Table);
  }
  /*   END: Names for output variables */
  for(i=0; i < SUB_OUTPUT_VARIABLES; i++)
    Table->OUTPUT_VARIABLE_INDEX[i] = Index_Output_Variables[i];
    // Up to i=(SUB_OUTPUT_VARIABLES-1) are set through command line values:
    // -n [SUB_OUTPUT_VARIABLES] 
  
  /* Some implementations of this code require to alloc memmory according to
     a number of state variables that can change dynamically. In that case, this
     3 lines of code should be commented out. Notice alse the corresponding lines
     in the function below:
     void P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( Parameter_Table * T )
     Only if the number of dynamic state model variables never cange during execution,
     these lines of code make sense here
  */
  /* Initial Conditions: MODEL STATE VARIABLES */
  // for (i=0; i < MODEL_STATE_VARIABLES; i++){
  //    AssignLabel_to_Model_Variables(i, Table->Model_Variable_Name[i], Table);
  // }

  assert(Table->TYPE_of_NETWORK == 1) ;  /* This ensures the network type is a squared grid 
                                            * with Von Neumann neighborhood 
                                            */
  Setting_up_Constant_Metapopulation_Connectivity_Matrix (Table);                                                                         
  
  /* This function should be called always after having called 
     void Parameter_Values_into_Parameter_Table(Parameter_Table * P)
  */
  #ifndef DIFFUSION_ECO_PLASMIDS
    Resetting_Lambda_Delta_Vectors (Table);
  #endif
  #ifndef DIFFUSION_ECO_1B1P
    Resetting_Lambda_Delta_Vectors (Table);
  #endif 
  /* END -------------------------------------------------*/
}

void Ressetting_Species_Characteristic_Parameters (Parameter_Table * Table)
{
  /* This function defines all species-specific parameters. 
     It is used to introduce tradeoffs between parameters 
  */

  int i; 
  
  assert(Table->TYPE_of_MODEL == 20);  /* 20: MODEL=DIFFUSION_ECOEVO_PLANTS */
  
  double p       = Table->p_1;
  double r_0     = Table->p_2;
  double z_0     = Table->p_2;
  double Eta_MIN = Table->Delta_C_0; 
  double Eta_MAX = Table->Delta_C_1;
   
  for(i=0; i<Table->No_of_RESOURCES; i++) {

    Table->Delta_RP[i] = Table->Delta_R_0; 
    Table->Delta_AP[i] = Table->Delta_R_1;
    
    Table->Eta_RP[i]= Eta_MIN + (double)i/((double)Table->No_of_RESOURCES - 1.0) * (Eta_MAX - Eta_MIN);
    /* Tradeoff: Between Eta and Beta based on an invasion criterion (r_0) */
    Table->Beta_AP[i] = r_0 * (Table->Eta_RP[i] + Table->Delta_RP[i])/(Table->Eta_RP[i]*(1.0-2.0*p)) * Table->Delta_AP[i]; 
    /* Tradeoff: Between Eta and Beta based on the amount of free space (z_0): equi-fitness manifold */
    // Table->Beta_AP[i]=(Table->Delta_RP[i] + z_0*Table->Eta_RP[i])/(z_0*Table->Eta_RP[i]) * Table->Delta_AP[i]; 

    Table->Mu_RP[i]  = Table->Mu; 
  }
}

void Resetting_Lambda_Delta_Vectors (Parameter_Table * Table)
{
  int i;

  if (Table->No_of_RESOURCES > 0 ) {
    Table->Lambda_R[0] = Table->Lambda_R_0;
    Table->Delta_R[0]  = Table->Delta_R_0;
  }
  if (Table->No_of_RESOURCES > 1 ) {
    Table->Lambda_R[1] = Table->Lambda_R_1;
    Table->Delta_R[1]  = Table->Delta_R_1;
    
  }
  if (Table->No_of_RESOURCES > 2 ) {
    for(i=2; i < Table->No_of_RESOURCES; i++) {
      Table->Lambda_R[i] = Table->Lambda_R_0;
      Table->Delta_R[i]  = Table->Delta_R_0;
    }
  }
}

void Resetting_Alpha_Nu_Vectors (Parameter_Table * Table)
{
  /* DIFFUSION_HII_nD: when a Holling Type II consumer feed on multiple resources */
  int i;

  if (Table->No_of_RESOURCES > 0 ) {
    Table->Alpha_C[0] = Table->Alpha_C_0;
    Table->Nu_Consumers[0]    = Table->Nu_C_0;
  }
  if (Table->No_of_RESOURCES > 1 ) {
    Table->Alpha_C[1] = Table->Alpha_C_0;  /* Overloading: */
    Table->Nu_Consumers[1]    = Table->Lambda_R_0; /* Nu_Consumers_1 would be Lambda_R_0 */
  }
  if (Table->No_of_RESOURCES > 2 ) {
    Table->Alpha_C[2] = Table->Alpha_C_0;  /* Overloading: */
    Table->Nu_Consumers[2] = Table->Lambda_R_1; /* Nu_Consumers_2 would be Lambda_R_1 */
  }
  if (Table->No_of_RESOURCES > 3 ) {
    for(i=3; i < Table->No_of_RESOURCES; i++) {
      Table->Alpha_C[i]  = Table->Alpha_C_0;
      Table->Nu_Consumers[i]  = Table->Nu_C_0;
    }
  }
}

void Resetting_Alpha_Nu_Vectors_Constant (Parameter_Table * Table)
{
  /* DIFFUSION_HII_nD: when a Holling Type II consumer feed on multiple resources */
  int i;
  for(i=0; i<Table->No_of_RESOURCES; i++) {
    Table->Nu_Consumers[i] = Table->Nu_C_0;
    Table->Alpha_C[i]      = Table->Alpha_C_0;
  }
}

void Resetting_Multiresource_Levels (Parameter_Table * Table)
{
  /* When a Holling Type II consumer feeds on multiple resources, and this multiple 
     resources are maintained, each of them at a constant density level
  */
  /* In this example, the different resource densities are arbitrarily fixed. 
     A total resource density is set to R (TOTAL_No_of_RESOURCES), and then:  
     1st Resource Density: 0.5 * R
     2on Resource Density: 0.3 * R
     For the rest of resources, they share the 0.2 * R that is left. 
  */
  int j;
  double K_R;

  K_R = (double)Table->K_R; 

  if ( Table->No_of_RESOURCES == 5) {
    Table->y_R_i[0]              = 0.9 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[0]    = Table->Alpha_C[0] * Table->y_R_i[0]/K_R;

    Table->y_R_i[1]              = 0.7 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[1]    = Table->Alpha_C[1] * Table->y_R_i[1]/K_R;

    Table->y_R_i[2]              = 0.5 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[2]    = Table->Alpha_C[2] * Table->y_R_i[2]/K_R;

    Table->y_R_i[3]              = 0.3 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[3]    = Table->Alpha_C[3] * Table->y_R_i[3]/K_R;

    Table->y_R_i[4]              = 0.1 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[4]    = Table->Alpha_C[4] * Table->y_R_i[4]/K_R;
  }
  else {
    if (Table->No_of_RESOURCES > 0 ) {
    Table->y_R_i[0]        = 0.5 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[0]    = Table->Alpha_C[0] * Table->y_R_i[0]/K_R;
    }
    if (Table->No_of_RESOURCES > 1 ) {
    Table->y_R_i[1]        = 0.3 * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[1]    = Table->Alpha_C[1] * Table->y_R_i[1]/K_R;
    }
    if (Table->No_of_RESOURCES > 2 ) {
    Table->y_R_i[2] = 0.2 / ((double)(Table->No_of_RESOURCES) - 2.0) * Table->TOTAL_No_of_RESOURCES;
    Table->Theta_Consumers[2] = Table->Alpha_C[2] * Table->y_R_i[2]/K_R;
    }
    if (Table->No_of_RESOURCES > 3 ) {
      for(j=3; j < Table->No_of_RESOURCES; j++) {
      Table->y_R_i[j]  = 0.2 / ((double)(Table->No_of_RESOURCES) - 2.0) * Table->TOTAL_No_of_RESOURCES;
      Table->Theta_Consumers[j]  = Table->Alpha_C[j] * Table->y_R_i[j]/K_R;
      }
    }
  }
}

void Writing_Alpha_Nu_Theta_Vectors(Parameter_Table * Table)
{
  int i; 
  double Density, Theta, T_C; 

  printf("Resource Type \t|\t Nu\t|\t Alpha\t|\t Theta\t|\t Density (y/K)\n");
  double Nu = Table->Nu_C_0; 
  
  Theta = 0.0; 
  for(i=0; i<Table->No_of_RESOURCES; i++) {
    
    Density = Table->y_R_i[i]/(double)Table->K_R;
    
    printf("%d:\t", i);
    printf("Nu_%d = %g\t", i, Table->Nu_Consumers[i]);
    printf("Alpha_%d = %g\t", i, Table->Alpha_C[i]);
    printf("Theta_%d = %g\t", i, Table->Theta_Consumers[i]);
    printf("Resource Density (fraction) y_R[%d]/K = %g\n", i, Density);

    Theta += Table->Theta_Consumers[i];
  }

  T_C = 1.0 / (Nu + Theta);

  Table->Tiempo_Propio = T_C; 
  
  printf(" Exact Characteristic Time (only is the handling time is the same across all resources):\n"); 
  printf(" T = %g\n", T_C);

  Print_Press_Key(0,0,".");
} 

/* void Parameter_Table_Index_Update(int * Index, int N, Parameter_Table * P) */
/* {                                                                          */
/*   int i;                                                                   */
/*   for(i=0; i<N; i++) P->Index[i] = Index[i];                               */
/* }                                                                          */
/* -------------------------------------------------------------------------- */
#if defined DIFFUSION_ECO_PLASMIDS
void Setting_Interaction_Matrices (Parameter_Table * Table) 
{
    int i, j;

    assert(Table->TYPE_of_MODEL == 21);  /* 21: MODEL=DIFFUSION_ECO_PLASMIDS */

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
  
  int n_bool_infection;
  int n_bool_incompatibility;  

  int No_of_POTENTIAL_PROFILES_per_STRAIN;
  int i, j, k, l; 

  No_of_POTENTIAL_PROFILES_per_STRAIN = pow(2.0, Table->No_of_PLASMIDS);

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

  N = 0;  
  for(i = 0; i < Table->No_of_STRAINS; i++) {

    Table->n_0[i] = N;  /* Saving the Strain ID of a plasmid free strain */

    n[i] = 0;  
    for(j = 0; j<No_of_POTENTIAL_PROFILES_per_STRAIN; j++) {
      
      n_bool_infection = 1; n_bool_incompatibility = 0;
      for(k = 0; k<Table->No_of_PLASMIDS; k++) {
        if( Profile[j][k] == 1) {  

        /* Infection Constraint */ 
          if(Table->IBP[i][k] == 1.0) n_bool_infection = 1; 
          else                        n_bool_infection = 0;
        
          n_bool_incompatibility = 0;
          for(l = 0; l<Table->No_of_PLASMIDS; l++) {
            if( Profile[j][l] == 1) {
            /* Plasmid Incompatibility Co-Infection Constraints */
              if(Table->CPP[k][l] == 0.0) 
                n_bool_incompatibility += 1;
            }
          } 
        }
      }

      if(n_bool_infection == 1 && n_bool_incompatibility == 0) { 

        for(k = 0; k<Table->No_of_PLASMIDS; k++)         
          Strain_Profiles[i][n[i]][k] = Profile[j][k];

        Table->StrainType_and_Profile[N][0] = i;     /* Strain Type          */
        Table->StrainType_and_Profile[N][1] = n[i];  /* Profile No Specifier */

        n[i]++;
        N++;                
      }      
    }        
  }

  printf("Number of different strains and profiles: %d \n", N);
  N = 0; 
  for(i = 0; i < Table->No_of_STRAINS; i++) {

    assert( Table->n_0[i] == N);
    
    N += n[i];
    
    for(j = 0; j < n[i]; j++) {
      printf("Strain Type [Profile] %d: ", i);

      printf("[ ");      
      for(k=0; k<Table->No_of_PLASMIDS; k++)
        printf("%d ", Strain_Profiles[i][j][k]);

      printf("]\n"); 
    }
    printf("\n");
  }  

  assert( N == Table->No_of_RESOURCES );

  for(i=0; i < No_of_POTENTIAL_PROFILES_per_STRAIN; i++) 
    free(Profile[i]);
  free(Profile);

  return(N); 
}

void Setting_Adjacency_Lists_from_Interaction_Matrices (Parameter_Table *Table)
{
  int i, j, k, N;
  int nc, nh, ni, nd;  
  int i_Focal, i_List;
  int k_Focal, k_List; 

  /* Adjacent Lists: */
  Table->No_of_CONJUGATION_EVENTS = 0;

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
    Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal);

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
    
    /* Recipient list of strain with the true index i (the DONOR)*/
    N = 0;
    for(k = 0; k<Table->No_of_PLASMIDS; k++)
     if(Table->Strain_Profiles[i_Focal][k_Focal][k] == 1)
      N++; 

    if(N == 0) { 
      assert( k_Focal == 0 );
      /* Empty Recipient List */
      Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES] = 0;
    }
    else {
      if( j != i ) {
        ni = 0; 
        for(j = 0; j<Table->No_of_RESOURCES; j++) {
          Calculate_Strain_and_Profile(Table, j, &i_List, &k_List);

          N = 0;
          if(Table->HBB[i_List][i_Focal] > 0.0)      /* Conjugation is permitted */ 
            for(k = 0; k<Table->No_of_PLASMIDS; k++) /* Infection is permitted */
              if(Table->Strain_Profiles[i_Focal][k_Focal][k] == 1 && Table->IBP[i_List][k] > 0.0) 
                N++;

          if(N > 0)
            Table->Recipient_List_Indeces[i][ni++] = j;
        }       
        Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES] = ni;
        Table->n_R[i]                                            = ni; /* No of Recipients of Strain ID 'i' */
        Table->No_of_CONJUGATION_EVENTS                         += ni;
      }         
    }
    
    /* Donor list of strain with the true index i (the RECIPIENT) */
    N = 0;
    for(k = 0; k<Table->No_of_PLASMIDS; k++)
     if(Table->Strain_Profiles[i_Focal][k_Focal][k] == 1)
      N++; 

    if(N == Table->No_of_PLASMIDS) { 
      /* Empty Donor List */
      Table->Donor_List_Indeces[i][Table->No_of_RESOURCES] = 0;
    }
    else {
      if( j != i ) {
        nd = 0; 
        for(j = 0; j<Table->No_of_RESOURCES; j++) {
          Calculate_Strain_and_Profile(Table, j, &i_List, &k_List);
          N = 0;
          if(Table->HBB[i_List][i_Focal] > 0.0)      /* Conjugation is permitted */ 
            for(k = 0; k<Table->No_of_PLASMIDS; k++) /* Infection is permitted   */
              if(Table->Strain_Profiles[i_List][k_List][k] == 1 && Table->IBP[i_Focal][k] > 0.0) 
                N++;
            
          if(N > 0)
            Table->Donor_List_Indeces[i][nd++] = j;      
        } 
        Table->Donor_List_Indeces[i][Table->No_of_RESOURCES] = nd;
      }  
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

void Calculate_Strain_and_Profile(Parameter_Table * Table, int n, int * i_Strain, int * k_Profile)
{
  /* Input: 
     . n, ID label of a subpopulation or Strain ID or local variable index 
          (true index of the local variable) 
     Output:
     . i_Strain, label of the strain type (from 0 to No_of_STRAINS-1)
     . k_Profile, label of the profile (within the i_Strain)
  */
  int i;
  int N; 

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
}

void Printing_Strains_Profiles_and_Lists(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 

  for(i=0; i<Table->No_of_STRAINS; i++)
    printf("Strain Type: %d\t No of Profiles of this strain type: %d\n", i, Table->n[i]);
  printf("\n");

  /* Adjacent Lists: */
  for(i = 0; i<Table->No_of_RESOURCES; i++) {
   Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 
   printf("Strain ID: %d\t Strain Type: %d\t Strain Profile [", i, i_Focal);
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

    printf("\n");
    printf("Recipient List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Recipient_List_Indeces[i][j]);
    printf("}"); 

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
#endif

#if defined(DIFFUSION_ECO_PLASMIDS) || defined(DIFFUSION_ECO_1B1P)
void Setting_Plasmid_Characteristic_Parameters (Parameter_Table * Table)
{
    //Plasmid Cost and Resistance 
    int i; 
    
    assert(Table->TYPE_of_MODEL == 21 || Table->TYPE_of_MODEL == 22);  /* 21: MODEL=DIFFUSION_ECO_PLASMIDS */
                                                                       /* 22: MODEL=DIFFUSION_ECO_1B1P    */
    for (i=0; i<Table->No_of_PLASMIDS; i++) {
      Table->Alpha_C[i]  = Table->Alpha_C_0;    /* Plasmid reproduction costs */
      Table->Nu_C[i]     = Table->Nu_C_0;       /* Plasmid resistance */
    }

} 

void Setting_Strain_Characteristic_Parameters (Parameter_Table * Table)
{
    int i, j, k, n;
    double COST, RESISTANCE;  
    double c;
    
    if(Table->TYPE_of_MODEL == 21) {  /* 21: MODEL=DIFFUSION_ECO_PLASMIDS */
  
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
    else if (Table->TYPE_of_MODEL == 22 ){             /* 22: MODEL=DIFFUSION_ECO_1B1P    */
      /* Plasmid Free Strain */
      Table->Beta_AP[0]  = Table->Beta_R;                                               /* Bacteria Cell Division Rate */
      Table->Eta_RP[0]   = Table->Lambda_R_1;                                           /* Bacteria Conjugation Rate   */
      Table->Delta_AP[0] = Table->Delta_R_0 + Table->Delta_R_1;                         /* Bacteria Death Rate         */
      Table->Mu_RP[0]    = Table->Mu;                                                   /* Bacteria Diffusion Rates    */
      Table->Segregation_Error[0]= 0.0;                                                 /* Bacterial Segregation Error */

      /* Plasmid Carrying Strain */
      Table->Beta_AP[1]  = Table->Beta_R * (1.0 - Table->Alpha_C[0]);                   /* Bacteria Cell Division Rate */
      Table->Eta_RP[1]   = Table->Lambda_R_1;                                           /* Bacteria Conjugation Rate   */
      Table->Delta_AP[1] = Table->Delta_R_0 + Table->Delta_R_1*(1.0 - Table->Nu_C[0]);  /* Bacteria Death Rate         */
      Table->Mu_RP[1]    = Table->Mu_C;                                                 /* Bacteria Diffusion Rates    */
      Table->Segregation_Error[1]= Table->p_1;                                          /* Bacterial Segregation Error */

      /* Notice that Alpha_C[0] is the cost of carrying the plasmid. 
         and Nu_C[0] is the resistance to the antibiotic. The pair (Alpha_C[0], Nu_C[0])
         characterizes the single plasmid in the system and defines a potential 
         cost-resistance tradeoff.  
      */
      assert(Table->No_of_RESOURCES == 2);      
    }
    else {
      printf(" This function can only be used with MODEL=DIFFUSION_ECO_PLASMIDS (TYPE_of_MODEL = 21)\n");
      printf(" or MODEL=DIFFUSION_ECO_1B1P (TYPE_of_MODEL == 22)\n");
      printf(" but the current model is TYPE_of_MODEL = %d\n", Table->TYPE_of_MODEL);
      assert(0);
    } 
}  
#endif

void Parameter_Values_into_Parameter_Table(Parameter_Table * P)
{
  /*
   The purpose of this simple function is just to upload
   ALL input parameters, controling both model definition
   and running simulations, which are defined as global variables,
   into the corresponding Parameter_Table Structure
  */ 
  P->Mu_C        = Mu_C; 
  P->Mu          = Mu;
  
  P->Lambda_R_0 = Lambda_R_0; 
  P->Delta_R_0  = Delta_R_0; 

  P->Lambda_R_1 = Lambda_R_1; 
  P->Delta_R_1  = Delta_R_1;

  P->K_R        = K_R;            /* -HK */
  P->Beta_R     = Beta_R;         /* -H4 */ 
  
  P->Lambda_C_0 = Lambda_C_0;     /* -H5 */
  P->Delta_C_0  = Delta_C_0;      /* -H6 */
  P->Lambda_C_1 = Lambda_C_1;     /* -H7 */
  P->Delta_C_1  = Delta_C_1;      /* -H8 */ 
      
  P->Alpha_C_0  = Alpha_C_0;      /* -H9 */
  P->Nu_C_0     = Nu_C_0;         /* -H10 */
  
  P->Chi_C_0    = Chi_C_0;        /* -H11 */
  P->Eta_C_0    = Eta_C_0;        /* -H12 */

  P->N_E        = N_E;            /* -H14 */ /* Number of Energy Levels */
  P->f          = f;              /* -H15 */ /* Fecundity: Number of Offspring Individuals */ 
  P->i_0        = i_0;            /* -H16 */ /* Energy Level at Maturity  */
  P->Beta_C     = Beta_C;         /* -H17 */ /* Consummer Reproduction Rate */
  P->k_E        = k_E;            /* -H18 */ /* 2* k_E is the resourse value in energy units */
  P->Theta_C    = Theta_C;        /* -H19 */ /* Energy loss rate for maintenance */
  P->p_1        = p_1;            /* -Hp1 */ /* Cooperation probability 1st position in the triplet */ 
  P->p_2        = p_2;            /* -Hp2 */ /* Cooperation probability 2on position in the triplet */ 

  P->Eta_R      = Eta_R;          /* -H20 */ /* Propagule establisment Rate */
  
  P->No_of_IC = No_of_IC;
  P->TYPE_of_INITIAL_CONDITION = TYPE_of_INITIAL_CONDITION;
  P->INITIAL_TOTAL_POPULATION  = INITIAL_TOTAL_POPULATION;

  P->RESCALING_INITIAL_TOTAL_POPULATION = RESCALING_INITIAL_TOTAL_POPULATION;

  P->TYPE_of_ERROR_FUNCTION = TYPE_of_ERROR_FUNCTION;
  P->No_of_ERROR_PARAMETERS = No_of_ERROR_PARAMETERS;

  P->Error_Parameter_0 = Error_Parameter_0;
  P->Error_Parameter_1 = Error_Parameter_1;

  P->Err_0 = Err_0;
  P->Err_1 = Err_1;
  P->Err_2 = Err_2;
  P->Err_3 = Err_3;
  P->Err_4 = Err_4;
  P->Err_5 = Err_5;
  P->Err_6 = Err_6;
  P->Err_7 = Err_7;
  P->Err_8 = Err_8;
  P->Err_9 = Err_9;
  P->Err_10 = Err_11;
  P->Err_11 = Err_11;
  P->Err_12 = Err_12;
  P->Err_13 = Err_13;
  P->Err_14 = Err_14;
  P->Err_15 = Err_15;
  P->Err_16 = Err_16;

  /* Definition of type of network */
  P->TYPE_of_NETWORK    = TYPE_of_NETWORK;
                               /* 0: Fully Connected Network
				                          1: Square Grid with Von Neuman neighborhood
				                          More network structures under construction 
			                          */
  P->No_of_NEIGHBORS  = No_of_NEIGHBORS;
  
#ifdef DIFFUSION_ECO_PLASMIDS
  
  P->No_of_STRAINS    = No_of_RESOURCES;   /* -HS 100, for instance */
  P->No_of_PLASMIDS   = No_of_INDIVIDUALS; /* -HN 10, for instance  */
  P->No_of_RESOURCES  = P->No_of_STRAINS * pow(2.0, P->No_of_PLASMIDS); /* Potential number of different species */
  P->No_of_PROFILES   = pow(2.0, P->No_of_PLASMIDS); /* Potential number of different profiles per strain */

  #elif defined DIFFUSION_ECO_1B1P

  P->No_of_STRAINS    = No_of_RESOURCES;   /* -HS 1, always, please */
  P->No_of_PLASMIDS   = No_of_INDIVIDUALS; /* -HN 1, always, please */
  P->No_of_RESOURCES  = P->No_of_STRAINS * pow(2.0, P->No_of_PLASMIDS); /* Potential number of different species */
  P->No_of_PROFILES   = pow(2.0, P->No_of_PLASMIDS); /* Potential number of different profiles per strain */

  assert(P->No_of_RESOURCES == 2 && P->No_of_PLASMIDS == 1);
#else
  
  P->No_of_RESOURCES    = No_of_RESOURCES; 

#endif
}
