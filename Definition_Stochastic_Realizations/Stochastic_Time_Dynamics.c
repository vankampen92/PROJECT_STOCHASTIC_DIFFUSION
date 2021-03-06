#include <MODEL.h>

/* Important Notes:
 *   . Table->Vector_Model_Int_Variables will store global variables across the patch system
 */

void S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S( int i,
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
      j_Good++; /* Counting good times */
      for(k=0; k < Table->SUB_OUTPUT_VARIABLES; k++){
	kk = Table->OUTPUT_VARIABLE_INDEX[k];
	value = definition_OutPut_Variables(kk, Table->Vector_Model_Variables,
					    Time->Time_Vector[0], Table);
	Time->summ[k][j] += value;
	Time->summ_var[k][j] += value * value;

#if defined CPGPLOT_REPRESENTATION
 	P->CPG->y_Time[k][j_Good] = value;
#endif
	Time->Variable[i][k][j]   = value;
	Table->Vector_Output_Variables[k] = value;
      }
#if defined CPGPLOT_REPRESENTATION
      if( FROZEN_SYSTEM == 0 ){ P->CPG->x_Time[j_Good]    = Time_Current; }
      else {            P->CPG->x_Time[j_Good]    = Time->Time_Vector[j]; }
#endif
      Time->time_DEF[j_Good] = Time_Current;
      Time->count[j]++;
      Time->Accumulated_Variable[i][j] = (double)new;

      Time_Initial = Time->Time_Vector[j-1];
      Time_Final   = Time->Time_Vector[j];

#if defined CPGPLOT_REPRESENTATION    /* Plotting Time evolution */
      /* BEGIN: Grafical Representation per SUCCESSFUL time step */
      C_P_G___S_U_B___P_L_O_T_T_I_N_G___n___P_L_O_T_S( Table->CPG->DEVICE_NUMBER,
						       1+i, j_Good, Table );

      /* GRID REPRESENTATION */
      Community_Scatter_Plot_Representation(Table, i, j);   
      // Press_Key(); 
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
    printf(" Total population across the system at current time (t = %g)\n", 
	   Time_Current );
    
    Print_Meta_Community_Patch_System (Table);
    Press_Key();
#endif

    
  }/* go further to the next time           */

  fclose(FP);
}
