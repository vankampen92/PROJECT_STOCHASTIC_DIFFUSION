#include <MODEL.h>

/// #defined SAVING_SLICES_TO_FILE

int generic_Function_Parameter_2Dim_Scan_Improved( Parameter_Table * P, 
						   int No_of_POINTS_1, int Input_Parameter_1,
						   int No_of_POINTS_2, int Input_Parameter_2,
						   double (* GENERIC_FUNCTION) (Parameter_Table *),
						   double * W_GRID, 
						   char * Scan_Output_File,
						   int X_LINEAR, int Y_LINEAR )
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

    if (Y_LINEAR == 0 ) {
      // Linear Scale
      Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min );
      Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX );
    }
    else if (Y_LINEAR == 1 ) {
      // Log Scale
      Value_0 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_min ));
      Value_1 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_2, S->Parameter_MAX ));
    }
    else {
      printf(" Y_LINEAR = %d, but it can only take 0/1 values\n", Y_LINEAR);
      printf(" The program will exit\n");
      Print_Press_Key(1,0,".");
      exit(0); 
    }
    
    Value = Value_0 + k * (Value_1 - Value_0)/(double)(No_of_POINTS_2 - 1);
    y_Data[k]= Value;

    if (Y_LINEAR == 0 ) 
      // Linear Scale
      AssignVectorEntry_to_Structure(P, Input_Parameter_2, y_Data[k] );
    else if ( Y_LINEAR == 1 ) 
      // Logr Scale
      AssignVectorEntry_to_Structure(P, Input_Parameter_2, pow(10.0, y_Data[k]) );
    else {
      printf(" Y_LINEAR = %d, but it can only take 0/1 values\n", Y_LINEAR);
      printf(" The program will exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }

    if( X_LINEAR == 0 ) {
      // Linear Scale
      Value_0 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min );
      Value_1 = Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX );
    }
    else if ( X_LINEAR == 1 ) {
      // Log Scale
      Value_0 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_min ));
      Value_1 = log10(Parameter_Model_into_Vector_Entry( Input_Parameter_1, S->Parameter_MAX ));
    }
    else {
      printf(" X_LINEAR = %d, but it can only take 0/1 values\n", Y_LINEAR);
      printf(" The program will exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
    }
    
    for( j = 0; j < No_of_POINTS_1; j++ ){
	
	Value = Value_0 + j * (Value_1 - Value_0)/(double)(No_of_POINTS_1 - 1);
	x_Data[j] = Value;
	
	if( X_LINEAR == 0 )
	  // Linear Scale:
	  AssignVectorEntry_to_Structure(P, Input_Parameter_1, x_Data[j] );
	else if ( X_LINEAR == 1 )
	  // Log Scale
	  AssignVectorEntry_to_Structure(P, Input_Parameter_1, pow(10.0, x_Data[j]) );
	else {
	  printf(" X_LINEAR = %d, but it can only take 0/1 values\n", Y_LINEAR);
	  printf(" The program will exit\n");
	  Print_Press_Key(1,0,".");
	  exit(0);
	}

	/* 
	   double pow(double x, double y);
	   This function returns the value of x raised to the power of y.
	*/
	
	// Most important code line in this function:

	// printf("k = %d\t j = %d ----------------------------------\t", k, j);  
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
  printf("\n From generic_Function_Parameter_Scan_Improved.c:\n  End of 2D scan successfully\n"); //getchar();

  /* BEGIN : Saving to File  */
  FILE * fp_0 = fopen ( Scan_Output_File, "w" );  
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

