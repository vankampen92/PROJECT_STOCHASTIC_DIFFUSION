#include <MODEL.h>

void Patch_System_Initialization (Community ** PATCH, Parameter_Table * Table, double * y_INI)
{
  /* This function is called from Initial_Conditions_Stochastic_Dynamics.c. This is the calling structure: 
  
    main() (main.c) 
          --->  M_O_D_E_L___S_T_O( &Table ); (MODEL_STO.c) 
                --->  S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S (...); (Stochastic_Time_Dynamics.c)
                      ---> Initial_Conditions_Stochastic_Dynamics(...); (Initial_Conditions_Stochastic_Dynamics.c) 
                           ---> Patch_System_Initialization (PATCH, Table, y_INI); (Patch_System_Initialization.c)
  */
  int i,j, k, S;
  int i_Strain, k_Profile; 
  double x, x_S;

  /* Populations are initialized in agreement with y_INI  */

  /* S is the number of variables required to define the state of a single patch */
  S = Table->LOCAL_STATE_VARIABLES;

  x = 0.0;
  for(i=0; i < S; i++) {

    x_S = 0; 
    for (j = 0; j < Table->No_of_CELLS; j++) {

      PATCH[j]->n[i] = (int)y_INI[i + j*S];

#if defined DIFFUSION_ECO_PLASMIDS
      PATCH[j]->Local_Strain_Population[i]->n  = PATCH[j]->n[i];
#endif 

      x_S += y_INI[i + j*S];
      x   += y_INI[i + j*S];

#if defined DIFFUSION_1R1C
      if (i == 0) printf("Initial Population[Species R] in patch %d: %d\n", j, PATCH[j]->n[i]);
      if (i == 1) printf("Initial Population[Species A] in patch %d: %d\n", j, PATCH[j]->n[i]);
      if (i == 2) printf("Initial Population[Species AR] in patch %d: %d\n", j, PATCH[j]->n[i]);
      if (i == 3) printf("Initial Population[Species ARA] in patch %d: %d\n", j, PATCH[j]->n[i]);
#endif
    }

#if defined DIFFUSION_1R1C
    if (i == 0) printf("Initial Population over all patches [Species R]: %g\n", x_S);
    if (i == 1) printf("Initial Population over all patches [Species A]: %g\n", x_S);
    if (i == 2) printf("Initial Population over all patches [Species RA]: %g\n", x_S);
    if (i == 3) printf("Initial Population over all patches [Species ARA]: %g\n", x_S);
#endif

    printf("Total Population Size at Time 0 (Species Type %d, across cells): %g\n", i, x_S);
  }

/* Certain models require the initialization of the quantity of empy space in every patch */
#if defined DIFFUSION_ECOEVO_PLANTS
  for (j = 0; j < Table->No_of_CELLS; j++) {
    x_S = 0;      
    for(i=0; i < S; i++) {
      x_S += PATCH[j]->n[i];
    }
    PATCH[j]->m_0 = (int)Table->K_R - (int)x_S; 
  }
#endif 

#if defined DIFFUSION_ECO_PLASMIDS
  for (j = 0; j < Table->No_of_CELLS; j++) {
    x_S = 0;      
    for(i=0; i < S; i++) {
      x_S += PATCH[j]->n[i];
    }
    PATCH[j]->m_0 = Table->K_R - x_S; 
  }

  for(i=0; i < Table->No_of_STRAINS; i++) {    
    for (j = 0; j < Table->No_of_CELLS; j++) {
      /* Local total number of individual cells per bacterial type (regardless profile) in PATCH 'j' */
      PATCH[j]->Bacterial_Type_Population[i]  = Bacterial_Type_Population_per_Cell (Table, i, j);
    }
  }
  for(i=0; i < Table->No_of_PLASMIDS; i++) {  
    for (j = 0; j < Table->No_of_CELLS; j++) {
      /* Local total number of individual bacterial cells carrying the same plasmid type in PATCH 'j' */ 
      PATCH[j]->Plasmid_Type_Population[i]  = Plasmid_Type_Population_per_Cell (Table, i, j);   

      PATCH[j]->Local_Plasmid_Population[i]->n = PATCH[j]->Plasmid_Type_Population[i];
    }
  }  
#endif 

#if defined DIFFUSION_HII_nD
  int TOTAL_No_of_HANDLING_CONSUMERS = 0;
  assert(Table->No_of_CELLS == 1);
  for(i=0; i < S; i++) {
    // printf("Initial Number of Handling Consumers (in the %d-Patch) on the %d-Resource Type: %d\n",
    //         0, i, PATCH[0]->n[i]);
    TOTAL_No_of_HANDLING_CONSUMERS += PATCH[0]->n[i];
  }
  Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0 = TOTAL_No_of_HANDLING_CONSUMERS;
  Table->TOTAL_No_of_HANDLING_CONSUMERS        = TOTAL_No_of_HANDLING_CONSUMERS;
  Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0     = Table->TOTAL_No_of_CONSUMERS - TOTAL_No_of_HANDLING_CONSUMERS;
  Table->TOTAL_No_of_FREE_CONSUMERS            = Table->TOTAL_No_of_CONSUMERS - TOTAL_No_of_HANDLING_CONSUMERS;
#else
  printf("\n");
  printf("Initial Total Population (for some Species Types at Time 0): %g\t", Table->INITIAL_TOTAL_POPULATION);
  pripntf("Total Community Size  (across Species): %g\n", x);
#endif

// printf("Patch system successfully initialized\n");

#if defined VERBOSE
    Print_Press_Key(1,1,"Patch system successfully initialized");
#endif
}