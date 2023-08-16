#include <MODEL.h>

/// #defined SAVING_SLICES_TO_FILE

int generic_Function_Parameter_2Dim_Scan( Parameter_Table * P, 
					  int No_of_POINTS_1, int Input_Parameter_1,
					  int No_of_POINTS_2, int Input_Parameter_2,
					  double (* GENERIC_FUNCTION) (Parameter_Table *),
					  double * W_GRID, 
					  char * generic_Function_Parameter_Scan_Output_File )
{
  double z_LOWER, z_INTER, z_UPPER;
  double Value, Value_0, Value_1;
  int i,j, k,n;
  Parameter_Space * S = P->S; 
  
  /* This function calculates a 2DIM scan of a GENERIC_FUNCTON which 
     depends on model parameters (Parameter_Table).
 
     This is done by scanning the parameter space defined by Input_Parameter_1 
     and Input_Parameter_2. The boundaries of the parameter domain under study 
     are defined in the corresponding  boundary_[TYPE_of_BOUNDARY].c file. 
     They can only be changed through changing that file and re-compiling again.

     Input parameters are labeled according to the input (model) parameters 
     labels as appear in all the assing_[].c functions. 

     The output of the probram generates the generic_Function_Parameter_Scan.dat 
     file, a three column file, (x, y, z), where z = FUNCTION (x, y) and a matrix 
     arranged in the array W_GRID[] 
  */

  /* BEGIN : Allocating memory for saving data to plot a bifurcation  * * * * * * */
  /*         diagram for each variable  * * * * * * * * * * * * * * * * * * * * * */  
  double      ** z_SOL  = (double **)malloc( No_of_POINTS_2 * sizeof(double *) );
  for( i = 0; i < No_of_POINTS_2; i++){
    z_SOL[i] = (double *)malloc( No_of_POINTS_1 * sizeof(double) );
  }
  double * x_Data  = (double *)malloc(No_of_POINTS_1 * sizeof(double) ); 
  double * y_Data  = (double *)malloc(No_of_POINTS_2 * sizeof(double) ); 
  /*   END : Allocating memory for saving dynamical data * * * * * */

  n = 0;
  for( k = 0; k < No_of_POINTS_2; k++ ) {

    // Linear Scale
    Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min );
    Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX );

    // Log Scale
    // Value_0 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min ));
    // Value_1 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX ));
    
    Value = Value_0 + k * (Value_1 - Value_0)/(double)(No_of_POINTS_2 - 1);
    
    y_Data[k]= Value; 
    // Log Scale
    // AssignVectorEntry_to_Structure(P, Input_Parameter_2, pow(10.0, y_Data[k]) );
    // Linear Scale
    AssignVectorEntry_to_Structure(P, Input_Parameter_2, y_Data[k] );

    // Linear Scale
    // Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min );
    // Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX );

    // Log Scale
    Value_0 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min ));
    Value_1 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX ));
  
    for( j = 0; j < No_of_POINTS_1; j++ ){
	
	Value = Value_0 + j * (Value_1 - Value_0)/(double)(No_of_POINTS_1 - 1);
	x_Data[j] = Value;
	// Log Scale:
	AssignVectorEntry_to_Structure(P, Input_Parameter_1, pow(10.0, x_Data[j]) );
	// Linear Scale:
	// AssignVectorEntry_to_Structure(P, Input_Parameter_1, x_Data[j] );

	// Most important code line in this function: 
	z_SOL[k][j]    = GENERIC_FUNCTION ( P );
	W_GRID[n++]    = z_SOL[k][j]; 

#if defined VERBOSE	
	printf(" n = %d\t x = %g\ty = %g\tz = %g\n", n, x_Data[j], y_Data[k], z_SOL[k][j] );
#endif 
      }

#if defined CPGPLOT_REPRESENTATION
      /* 2-DIM bifurcation diagram */
      /* C_P_G___P_L_O_T_T_I_N_G___S_C_A_N ( P, No_of_POINTS_1,        */
      /* 					  x_Data, z_SOL[k],    */
      /* 					  Input_Parameter_1,   */
      /* 					  Input_Parameter_2 ); */
      //Print_Press_Key(1,0,".");
#endif
#if defined SAVING_SLICES_TO_FILE 
      Saving_to_File_double("Parameter_Scan_Slice_", x_Data, z_SOL[k], No_of_POINTS_1, k); 
#endif
  }
  printf("\n From generic_Function_Parameter_Scan.c:\n  End of 2D scan successfully"); getchar();

  /* BEGIN : Saving to File  */
  FILE * fp_0 = fopen ( generic_Function_Parameter_Scan_Output_File, "w" );  
  for( k = 0; k < No_of_POINTS_2; k++ ) {
      for( j = 0; j < No_of_POINTS_1; j++ ){
	fprintf(fp_0, "%g\t%g\t%g\n", x_Data[j], y_Data[k], z_SOL[k][j]);
      }
  }
  fclose(fp_0); ;
  /*   END : Saving to File */  
  
  /* BEGIN : Freeing previous allocated memory */
  for(i = 0; i < No_of_POINTS_2; i++) { 
    free (z_SOL[i]);
  }
  free(z_SOL); 

  free(x_Data); free(y_Data);
  /*   END : End freeing allocated memory */

  return(0);
}

#if defined CPGPLOT_REPRESENTATION
void C_P_G___P_L_O_T_T_I_N_G___S_C_A_N ( Parameter_Table * P, 
					 int NO_of_POINTS, double * x_Data, double * y_SOL, 
					 int Input_Parameter_1, 
					 int Input_Parameter_2 )
{
  /* 
     CPG Representation: This function represents a bifurcation diagram for Output Variable, 
     where the bifurcation parameter, represented in the x axis, is given as Input_Parameter. 
  */

  int i, SCALE_FIXED;
  char * p_Title;
  char * Y_label  = (char *)malloc( sizeof(char) * 500 );
  char * X_label; 
  char * Title_In = (char *)malloc( sizeof(char) * 500 );
  char * Number   = (char *)malloc( sizeof(char) * 50  );
  /*********************************************************************/
  X_label   = P->Name_Parameters[Input_Parameter_1]; 
  Y_label[0] = '\0';
  p_Title = strcat( Y_label, "  z = F ( x ) ");
  /*********************************************************************/
  
  Title_In[0] = '\0';
  p_Title = strcat( Title_In, "  P a r a m e t e r   S c a n  ");
  p_Title = strcat( Title_In, "  ( ");
  p_Title = strcat( Title_In, P->Name_Parameters[Input_Parameter_2] );
  p_Title = strcat( Title_In, "  =  ");
  
  double Value = AssignStructValue_to_VectorEntry( Input_Parameter_2, P );
  doubletochar( Value, Number );

  p_Title = strcat( Title_In, Number);
  p_Title = strcat( Title_In, " )  ");
  
  //SCALE_FIXED  = 0;
  int SCALE_X = 0; 
  int SCALE_Y = 0;
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_C_A_L_E ( P->CPG,
						NO_of_POINTS, x_Data, y_SOL,
						X_label, Y_label, Title_In,
						SCALE_X, SCALE_Y );

  // C_P_G___2d___X_Y___P_L_O_T_T_I_N_G ( NO_of_POINTS, x_Data, y_SOL, 
  //				       X_label, Y_label, Title_In );
  
  free (Title_In); free (Y_label); free (Number);
} 

void Parameter_Scan_2D_GRID_SHADES( Parameter_Table * Table,
				    int No_of_POINTS_1, int Input_Parameter_1, 
				    int No_of_POINTS_2, int Input_Parameter_2, 
				    double (* GENERIC_FUNCTION) (Parameter_Table *), 
				    char * Output_File_Data )

  {  
    int i;
    double Value_0, Value_1;
    double * W_GRID = (double *)malloc( No_of_POINTS_1 * No_of_POINTS_2 * sizeof(double) );

    int Status =  generic_Function_Parameter_2Dim_Scan( Table,
							No_of_POINTS_1, Input_Parameter_1,
							No_of_POINTS_2, Input_Parameter_2,
							GENERIC_FUNCTION,
							W_GRID,
							Output_File_Data );

    Parameter_Space * Space = Table->S;
    
    /* BEGIN : 2D GRID cpgplot representation */
    /*********************************************************************/
    Table->CPG->X_label   = Table->Name_Parameters[Input_Parameter_1]; 
    Table->CPG->Y_label   = Table->Name_Parameters[Input_Parameter_2]; 
    /*********************************************************************/
    // Boundary(Input_Parameter_1, &Value_0, &Value_1);
    Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_min );
    Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, Space->Parameter_MAX );
    
    Table->CPG->ORIGIN_X    = Value_0;
    Table->CPG->X_Dimension = (Value_1 - Value_0);
    
    // Boundary(Input_Parameter_2, &Value_0, &Value_1);
    Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_min );
    Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, Space->Parameter_MAX );
    
    Table->CPG->ORIGIN_Y = Value_0;
    Table->CPG->Y_Dimension = (Value_1 - Value_0);
    
    Table->CPG->x_GRID  = No_of_POINTS_1; 
    Table->CPG->y_GRID  = No_of_POINTS_2;
    
    // int Output_Variable  = Table->IO_VARIABLE_LIST[0];
    // Table->CPG->W_label   = Table->Variable_Name[Output_Variable];
    
    int FIRST_PLOT = 0;
    double i_PLOT  = 0.0;
    
    // C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___S_H_A_D_E_S ( ... );
    // C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___C_O_N_T_O_U_R( ... ); 
    
    double W_GRID_MAX = 0.0;
    for(i=0; i<(No_of_POINTS_1 * No_of_POINTS_2); i++){
      W_GRID_MAX = MAX ( W_GRID[i], W_GRID_MAX );
    }  
    for(i=0; i<(No_of_POINTS_1 * No_of_POINTS_2); i++){
      if(W_GRID[i] > 20.0) W_GRID[i] = 20.0;
    }  
    
    C_P_G___P_L_O_T_T_I_N_G___2d___G_R_I_D___S_H_A_D_E_S( Table->CPG,
							  W_GRID, 
							  FIRST_PLOT,
							  Table->CPG->CPG_SCALE_W, 
							  Table->CPG->CPG_RANGE_W_0,
							  Table->CPG->CPG_RANGE_W_1,
							  i_PLOT );
    free (W_GRID);
  }
#endif
