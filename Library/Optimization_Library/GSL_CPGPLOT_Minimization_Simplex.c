#include <MODEL.h>



void GSL_CPGPLOT_Minimization_Simplex (Parameter_Fitting * F,
				       gsl_vector * Solution, size_t iter,
				       double ( * Function )( const gsl_vector * , void * ) )
{
  #if defined CPGPLOT_REPRESENTATION
  int i, key;
  double value;
  static int Iteration;

  if (iter == 0) Iteration++;
  else Iteration = (int)iter;
  
  Parameter_Table * Table  = F->Table;
  Parameter_Space * Space  = F->Space;

  value = ( * Function )(Solution, F);

  if (value == DBL_MAX )
    printf(" GSL_CPGPLOT_Minimization_Driver: Function evaluation outside boundaries\n");
  else {
    /* B E G I N : Creating nice title with parameter values and greek symbols */
    double Parameter_Value;
    char * pTitle;
    char * pValue = (char *)calloc(20, sizeof(char) );
    char * Title_Parameters = (char *)calloc(2000, sizeof(char) );
    Title_Parameters[0]='\0';
    pTitle = strcat(Title_Parameters, F->Data->Name);
    pTitle = strcat(Title_Parameters, " Iter: ");
    pValue[0]='\0'; sprintf(pValue, "%d", Iteration);
    pTitle = strcat(Title_Parameters, pValue);
    pTitle = strcat(Title_Parameters, ": ");
    for(i=0; i<Space->No_of_PARAMETERS; i++) {
      key = Space->Parameter_Index[i];
      pTitle = strcat(Title_Parameters, Table->Symbol_Parameters_Greek_CPGPLOT[key]);
      pTitle = strcat(Title_Parameters, " = ");
      Parameter_Value  = AssignStructValue_to_VectorEntry ( key, Table );
      pValue[0]='\0'; sprintf(pValue, "%4.3f", Parameter_Value);
      pTitle = strcat(Title_Parameters, pValue);
      pTitle = strcat(Title_Parameters, "  ");  
    }
    for(i=0; i<Table->No_of_IC; i++) {
      key = Table->IC_Space->Parameter_Index[i];
      pTitle = strcat(Title_Parameters, Table->Model_Variable_Symbol[key]);
      pTitle = strcat(Title_Parameters, " = ");
      Parameter_Value  = Model_Variable_Initial_Condition_into_Vector_Entry_Table( key, Table );
      pValue[0]='\0'; sprintf(pValue, "%4.3f", Parameter_Value);
      pTitle = strcat(Title_Parameters, pValue);
      pTitle = strcat(Title_Parameters, "  ");  
    }
    for(i=0; i<Table->No_of_ERROR_PARAMETERS; i++ ) {
      key = Table->E_Space->Parameter_Index[i];
      pTitle = strcat(Title_Parameters, "Err(");
      if(key < Table->OUTPUT_VARIABLES_GENUINE)
	pTitle = strcat(Title_Parameters, Table->Output_Variable_Symbol[key]);
      pTitle = strcat(Title_Parameters, ") = ");
      Parameter_Value  = Error_Model_into_Vector_Entry_Table( key, Table );
      pValue[0]='\0'; sprintf(pValue, "%4.3f", Parameter_Value);
      pTitle = strcat(Title_Parameters, pValue);
      pTitle = strcat(Title_Parameters, "  ");     
    }
    pTitle = strcat(Title_Parameters, ". NegLogLK:  ");
    if (value < 1.E+18) {  
      pValue[0]='\0'; sprintf(pValue, "%4.3f", value);
      pTitle = strcat(Title_Parameters, pValue);
    }
    else {
      pTitle = strcat(Title_Parameters, " inf");
    }
    pTitle = strcat(Title_Parameters, ".");
    /*     E N D : ----------------------------------------- */

    printf("Title of Plot: %s\n", Title_Parameters); 
    
    int SAME, k;
    for(k = 0; k<Table->SUB_OUTPUT_VARIABLES; k++) {

      key  = Table->OUTPUT_VARIABLE_INDEX[k];

      Table->CPG->CPG_SCALE_X   = 1;        Table->CPG->CPG_SCALE_Y   =  1;
      Table->CPG->CPG_RANGE_X_0 = Table->T->Time_Vector[0];
      Table->CPG->CPG_RANGE_X_1 = Table->T->Time_Vector[Table->T->I_Time-1];

      SAME = 0;
      Table->CPG->color_Index = 5; //* 5: green *//
      Table->CPG->type_of_Line = 1;
      Table->CPG->type_of_Width = 7;
      Table->CPG->type_of_Symbol = 1;

      CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( Table->CPG, SAME,
							    Table->T->I_Time,
							    Table->T->Time_Vector,
							    Table->Matrix_Output_Variables[k],
							    "Time",
							    Table->Output_Variable_Name[key],
							    Title_Parameters,
							    Table->CPG->CPG_SCALE_X,
							    Table->CPG->CPG_SCALE_Y );

      SAME = 1;
      Table->CPG->color_Index = 8; //* 3: green *//
      Table->CPG->type_of_Line = 1;
      Table->CPG->type_of_Width = 7;
      Table->CPG->type_of_Symbol = 24;

      if(k < F->Data->No_of_VARIABLES)
      CPGPLOT___X_Y___S_C_A_T_T_E_R_E_D___S_A_M_E___P_L_O_T ( Table->CPG, SAME,
							      F->Data->No_of_POINTS, 
							      Table->T->Time_Vector,
							      F->Data->N[k],
							      "Year",
							      Table->Output_Variable_Name[key],
							      Title_Parameters,
							      Table->CPG->CPG_SCALE_X,
							      Table->CPG->CPG_SCALE_Y );
    }
    free (pValue);
    free(Title_Parameters);
  }
  #if defined VERBOSE
  // Print_Press_Key(1,0,".");
  #endif
  
  #endif
  
}
