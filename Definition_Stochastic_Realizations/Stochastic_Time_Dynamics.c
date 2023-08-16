#include <MODEL.h>

/* Important Notes:
 *   . Table->Vector_Model_Int_Variables will store global variables across the patch system
 */

// #define EXTINCTION_CONTROL

int S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S( int i,
		  				                                       Parameter_Table * Table,
						                                         int * Bad_Times )
{
  /* This function perform one single stochastic realization (i-th),
     sample the system at times stored in Time->Time_Vector[],
     and save a file corresponding to this i-th stochastic
     realization in re_[i].dat file.
  */
  FILE *FP; char file[12];
  int j, k, kk, j_Good, Sp;
  int new; /* Ever-increasing Accumulated Variable within a time interval */
  int TIMES;
  Time_Control * Time;
  Parameter_Model * P;
  Community ** PATCH;
  double Time_Initial, Time_Current, Time_Final, value;
  Stochastic_Rate * Rate;

  P            = Table->P;
  PATCH        = Table->Patch_System;

  Time         = Table->T;
  Time_Initial = Time->Time_0;
  Time_Final   = Time->Time_1;
  TIMES        = Time->I_Time;
  Rate         = Time->Rate;

  /* Each stochastic realization will be saved in a different file */
  file[0]='\0';  fitxer(file, "re_", i, ".dat"); FP = fopen(file, "w");

  /* BEGIN : Initial Conditions -------------------------------------------------------------*/
  
  // printf(" Before  Initial_Conditions_Stochastic_Dynamics(...)\n");
  Initial_Conditions_Stochastic_Dynamics( Table, Table->Vector_Model_Variables );
  // printf(" After Initial_Conditions_Numerical_Integration(...). Initial Conditions:  ");

  Time_Current = Time->Time_Vector[0];
  Time->Time_Vector_Real[i][0] = Time_Current;  

  if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) {
    /* Update_Time_Dependence (Table->TYPE_of_TIME_DEPENDENCE, Time_Current, Table ); */
    // Time_Dependence_Apply( Table, Time_Current );
  }

#if defined CPGPLOT_REPRESENTATION
  P->CPG->x_Time[0]      = Time->Time_Vector[0];
#endif

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    kk = Table->OUTPUT_VARIABLE_INDEX[k];
    value = definition_OutPut_Variables(kk,
					                              Table->Vector_Model_Variables,
					                              Time->Time_Vector[0], Table);
#if defined CPGPLOT_REPRESENTATION
    P->CPG->y_Time[k][0] = value;
#endif
    Time->Variable[i][k][0]           = value;
    Table->Vector_Output_Variables[k] = value;
  }
  /* Initial calculation of the system total rate of a configurational change and rates of
     the different configurational changes or events to occur
  */
  Temporal_Dynamics(PATCH, Table, Rate);

  #if defined BINARY_TREE_OPTIMIZATION
    /* Initiate values of the binary tree with total rates of every patch at the leaves 
     */
    Community_Binary_Tree_Initialization (Table);
    P->Treeroot = Table->Treeroot;
    printf(" Binary Tree (from Leaves) to Sample Discrte Distribution has been successcully\n");
    printf(" initiated [Realization: %d]\n", i ); 
  #endif         
  /*   END : Initial Conditions -------------------------------------------------------------*/

  /* int DISCARTING_EXTINCTIONS = P->DISCARTING_EXTINCTIONS;  */
  int FROZEN_SYSTEM               = 0;
  (*Bad_Times)                    = 0; j_Good = 0;
  for( j = 1; j < TIMES; j++ ) {
    /* This loop advances the system sequentially from
       intitial time 0 to 1st time , ...,  from time (j-1) to j, and so on.
       Note: When the system is frozen (FROZEN_SYSTEM = 1), then
             this loop does not advance the system any more
    */
    /* B E G I N :
     *     CENTRAL POINT HERE: Stochastic Dynamics Loop While (up to the next time)
     */
    new = 0;
    while( Time_Current < Time->Time_Vector[j] && FROZEN_SYSTEM == 0 )
      {
	      FROZEN_SYSTEM = Advance_Current_Time( Table, Rate, &Time_Current, &new );
      }
    /*     E N D
     * -------------------------------------------------------------------
     */    
#if defined EXTINCTION_CONTROL
    int EXTINCTION;
    double Total_Population_of_Consumers;  
    EXTINCTION = 0;

    if (Table->TYPE_of_MODEL == 2 || Table->TYPE_of_MODEL == 8) {
      Total_Population_of_Consumers = Total_Population_Consumers (Table->Vector_Model_Variables, Table );
      if (Total_Population_of_Consumers == 0.0) EXTINCTION = 1;
    
    // EXTINCTION = Extinction_Control_Condition(Table);
    }
    
    if (FROZEN_SYSTEM == 1 || EXTINCTION == 1) {
      if (FROZEN_SYSTEM == 1 && EXTINCTION == 0) FROZEN_SYSTEM = 1;
      if (FROZEN_SYSTEM == 0 && EXTINCTION == 1) FROZEN_SYSTEM = 2;
      if (FROZEN_SYSTEM == 1 && EXTINCTION == 1) FROZEN_SYSTEM = 1;
      break;  /* Don't advance more times */  
    }
#endif

    if( Time_Current > Time->Time_Vector[j] + Time->EPSILON){
      (*Bad_Times)++;
#if defined VERBOSE
      printf("Time too far away from target: skipping this time...\n");
      printf("%d\t%g\t%g\n",j, Time_Current, Time->Time_Vector[j]);
#endif
    }
    else{
      /* Saving and representing at values close to Time_values[j] * * * * * * * * * * * * */
      /* Notice that Time_Current is always the last time which is the closest possible to
	      (and a little bit larger than) the time stored in Time->Time_Vector[j].           */
      j_Good++;            /* Counting good times of the i-th realization  */
      Time->count[j]++;    /* Counting good realizations corresponding to the j-th time */ 

      for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
	      kk = Table->OUTPUT_VARIABLE_INDEX[k];
	      value = definition_OutPut_Variables(kk, Table->Vector_Model_Variables, Time->Time_Vector[0], Table);
	      Time->summ[k][j] += value;
	      Time->summ_var[k][j] += value * value;

        #if defined CPGPLOT_REPRESENTATION
 	        P->CPG->y_Time[k][j_Good] = value;
        #endif
	      Time->Variable[i][k][j]   = value;
	      Table->Vector_Output_Variables[k] = value;
      }
      Time->Time_Vector_Real[i][j] = Time_Current; 
      
      #if defined CPGPLOT_REPRESENTATION
        if( FROZEN_SYSTEM == 0 ){ P->CPG->x_Time[j_Good]    = Time_Current; }
        else {            P->CPG->x_Time[j_Good]    = Time->Time_Vector[j]; }
      #endif
      
      Time->time_DEF[j_Good] = Time_Current;
      
      Time->Accumulated_Variable[i][j] = (double)new;

      Time_Initial = Time->Time_Vector[j-1];
      Time_Final   = Time->Time_Vector[j];

#if defined CPGPLOT_REPRESENTATION    /* Plotting Time evolution */
      /* BEGIN: Grafical Representation per SUCCESSFUL time step */
      C_P_G___S_U_B___P_L_O_T_T_I_N_G___n___P_L_O_T_S(Table->CPG->DEVICE_NUMBER,
      						                                    1+i, j_Good, Table );
      if( Table->No_of_CELLS > 4 ) 
        /* GRID REPRESENTATION */
	      Community_Scatter_Plot_Representation(Table, i, j);   
        // Print_Press_Key(1,0,"."); 
      /*   END: Grafical Representation per time step */
#endif
      /* BEGIN : Writing a costumized file ... */
      fprintf(FP,"%g", Time_Current);
      for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
	      fprintf(FP,"\t%g", Table->Vector_Output_Variables[k]);
      }
      fprintf(FP,"\n");
      /*   END: Writing costumized file        */
    }

#if defined VERBOSE
    printf(" Total population across the system at current time (t = %g)\n", Time_Current );
    Print_Meta_Community_Patch_System (Table);
    Print_Press_Key(1,0,".");
#endif
  }/* go further to the next time           */

#if defined BINARY_TREE_OPTIMIZATION
  // deleteBinaryTree_DiscreteDistribution(Table->Treeroot); 
  // Delete the Tree from root but not Leaves!!! 
#endif

  fclose(FP);

  return(FROZEN_SYSTEM); 
}

int Stochastic_Time_Dynamics_Numerical( int i,
					Parameter_Table * Table,
					int * Bad_Times )
{
  /* This version of the same function does the same as before, this is, 
     it performs one single stochastic realization (the i-th one). 
     No saving into files nor visually representing anything.  
  */
  int j, k, kk, j_Good, Sp;
  int new; /* Ever-increasing Accumulated Variable within a time interval */
  int TIMES;
  Time_Control * Time;
  Parameter_Model * P;
  Community ** PATCH;
  double Time_Initial, Time_Current, Time_Final, value;
  Stochastic_Rate * Rate;

  P            = Table->P;
  PATCH        = Table->Patch_System;

  Time         = Table->T;
  Time_Initial = Time->Time_0;
  Time_Final   = Time->Time_1;
  TIMES        = Time->I_Time;
  Rate         = Time->Rate;

    /* BEGIN : Initial Conditions -------------------------------------------------------------*/
  // printf(" Before  Initial_Conditions_Stochastic_Dynamics(...)\n");
  Initial_Conditions_Stochastic_Dynamics( Table, Table->Vector_Model_Variables );
  // printf(" After Initial_Conditions_Numerical_Integration(...). Initial Conditions:  ");

  Time_Current = Time->Time_Vector[0];
  Time->Time_Vector_Real[i][0] = Time_Current; 

  if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) {
    /* Update_Time_Dependence (Table->TYPE_of_TIME_DEPENDENCE, Time_Current, Table ); */
    // Time_Dependence_Apply( Table, Time_Current );
  }

  for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
    kk = Table->OUTPUT_VARIABLE_INDEX[k];
    value = definition_OutPut_Variables(kk,
					                              Table->Vector_Model_Variables,
					                              Time->Time_Vector[0], Table);
    Time->Variable[i][k][0]           = value;
    Table->Vector_Output_Variables[k] = value;
  }
  /* Initial calculation of the system total rate of a configurational change and rates of
     the different configurational changes or events to occur
  */
  Temporal_Dynamics(PATCH, Table, Rate);
  /* Set upf of the intial values of the binary tree with the total rates for every patch at 
     the leaves 
  */
  #if defined BINARY_TREE_OPTIMIZATION
    Community_Binary_Tree_Initialization (Table);
    P->Treeroot = Table->Treeroot;
  #endif
  /*   END : Initial Conditions -------------------------------------------------------------*/

  /* int DISCARTING_EXTINCTIONS = P->DISCARTING_EXTINCTIONS;   */
  int FROZEN_SYSTEM               = 0;
  (*Bad_Times)                    = 0; j_Good = 0;
  
  for( j = 1; j < TIMES; j++ ) {
    /* This loop advances the system sequentially from
       intitial time 0 to 1st time , ...,  from time (j-1) to j, and so on.
       Note: When the system is frozen (FROZEN_SYSTEM = 1), then
             this loop does not advance the system any more
    */
    /* B E G I N :
     *     CENTRAL PROGRAM CORE HERE: Stochastic Dynamics Loop While (up to the next time)
     */
    new = 0;
    while( Time_Current < Time->Time_Vector[j] && FROZEN_SYSTEM == 0 )
      {
	      FROZEN_SYSTEM = Advance_Current_Time( Table, Rate, &Time_Current, &new );
      }
    /*     E N D
     * -------------------------------------------------------------------
     */    
#if defined EXTINCTION_CONTROL
    int EXTINCTION;
    double Total_Population_of_Consumers;  
    EXTINCTION = 0;

    if (Table->TYPE_of_MODEL == 2 || Table->TYPE_of_MODEL == 8) {
      Total_Population_of_Consumers = Total_Population_Consumers (Table->Vector_Model_Variables,
								  Table );
      if (Total_Population_of_Consumers == 0.0) EXTINCTION = 1;
    
    // EXTINCTION = Extinction_Control_Condition(Table);
    }
    
    if (FROZEN_SYSTEM == 1 || EXTINCTION == 1) {
      if (FROZEN_SYSTEM == 1 && EXTINCTION == 0) FROZEN_SYSTEM = 1;
      if (FROZEN_SYSTEM == 0 && EXTINCTION == 1) FROZEN_SYSTEM = 2;
      if (FROZEN_SYSTEM == 1 && EXTINCTION == 1) FROZEN_SYSTEM = 1;
      break;  /* Don't advance more times */  
    }
#endif

    if( Time_Current > Time->Time_Vector[j] + Time->EPSILON){
      (*Bad_Times)++;
#if defined VERBOSE
      printf("Time too far away from target: skipping this time...\n");
      printf("%d\t%g\t%g\n",j, Time_Current, Time->Time_Vector[j]);
#endif
    }
    else{ 
      j_Good++; /* Counting good times */
      for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
	      kk = Table->OUTPUT_VARIABLE_INDEX[k];
	      value = definition_OutPut_Variables(kk, Table->Vector_Model_Variables, Time->Time_Vector[0], Table);
	      Time->summ[k][j] += value;
	      Time->summ_var[k][j] += value * value;

	      Time->Variable[i][k][j]   = value;
	      Table->Vector_Output_Variables[k] = value;
      }
      Time->Time_Vector_Real[i][j] = Time_Current; 

      Time->time_DEF[j_Good] = Time_Current;
      Time->count[j]++;
      Time->Accumulated_Variable[i][j] = (double)new;

      Time_Initial = Time->Time_Vector[j-1];
      Time_Final   = Time->Time_Vector[j];
    }
#if defined VERBOSE
    printf(" Total population across the system at current time (t = %g)\n",  Time_Current );
    Print_Meta_Community_Patch_System (Table);
    Print_Press_Key(1,0,".");
#endif
  }/* go further to the next time           */

  return(FROZEN_SYSTEM); 
}
