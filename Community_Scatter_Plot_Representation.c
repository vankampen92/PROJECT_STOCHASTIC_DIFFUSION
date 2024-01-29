#include <MODEL.h>
#include "Community.h"

extern gsl_rng * r; /* Global generator defined in main.c */

void Community_Scatter_Plot_Representation_4Sp( Parameter_Table * Table,
					    int i_Replicate, int j_Time )
{
  int i,j,k;
  int N, n;
  double x, x_0, c_x, y, y_0, c_y;
  static int DEVICE_NUMBER;
  Parameter_CPGPLOT * CPG;
  int Sp;
  float Range_x[2], Range_y[2];
  int color_Index;
  int type_of_Line;
  int type_of_Width;
  int type_of_Symbol;

  Community * Patch; 
  Community ** P = Table->Patch_System;
  double * Y = Table->Vector_Model_Variables;

  N =  Total_Population(Y, Table);

  if (Table->TYPE_of_MODEL == 0 )
    assert( N == (Table->No_of_RESOURCES * Table->No_of_INDIVIDUALS) );

  float * x_Data = (float *)calloc( N, sizeof(float) );
  float * y_Data = (float *)calloc( N, sizeof(float) );

  Range_x[0] = 0.0;  Range_x[1] = P[0]->X_DIMENSION;
  Range_y[0] = 0.0;  Range_y[1] = P[0]->Y_DIMENSION;

  if(i_Replicate == 0 && j_Time == 0) {
    DEVICE_NUMBER = cpgopen( "/TPNG" ); // "/XSERVE"
    cpgsubp(1, 1);
    cpgask( 0 );
  }

  for(Sp = 0; Sp<Table->LOCAL_STATE_VARIABLES; Sp++) {

    n=0;
    for(i = 0; i<Table->No_of_CELLS; i++) {
      Patch = P[i];
      
      for(j = 0; j<Patch->n[Sp]; j++) {
          x_Data[n] = (float)(Patch->center.x - 0.5 + gsl_rng_uniform(r));
	        /* This is because STEP is 1.0 */
	        y_Data[n] = (float)(Patch->center.y - 0.5 + gsl_rng_uniform(r));
	        /* This is because STEP is 1.0 */
	        n++;
      }
    }
      
    type_of_Line = 1;
    type_of_Width = 1;
    cpgsls(type_of_Line);
    cpgslw(type_of_Width);
    cpgask( 0 );

    cpgslct(DEVICE_NUMBER);

    if (Sp == 0) {
	    color_Index    = 2;        /* 2: Red        */ 
	    type_of_Symbol = 1;        /* 1: .  (point) */
    
      cpg_XY_scattered(n, x_Data, y_Data, Range_x, Range_y,
		       color_Index, type_of_Symbol,
		       "X", "Y", "");
    }
    else if (Sp == 1){
      color_Index    = 7;          /* 7: Yellow     */
	    type_of_Symbol = 2;          /* 2:  +         */
	    cpg_XY_same_scattered(n, x_Data, y_Data, 
			                      color_Index, type_of_Symbol);
    }
    else {
      color_Index    = 3 + Sp;    /* 4: Blue  5: Light Blue  6: Magenta    */
	    type_of_Symbol = 1 + Sp;    /* 2: +     3: *           4: o (cercle) */ 
	    cpg_XY_same_scattered(n, x_Data, y_Data, 
			                      color_Index, type_of_Symbol);
    }
    cpgsls(1);
    cpgslw(1);
  }
  
  /* n: Total number of individuals across species */
  // Print_Meta_Community_Patch_System (Table);

  /* B E G I N : Plotting Current Time on Top */
  char * Plot_Time  = (char *)calloc( 50, sizeof(char));
  char * Time_Eraser = (char *)calloc(50, sizeof(char));
  float char_Size;
  float x_Time_Position = 0.50 * Range_x[1]; 
  float y_Time_Position = 1.09 * Range_y[1];
  static double Current_Time  = 0.0; 
  double Last_Time            = Current_Time;
  Current_Time  = Table->T->time_DEF[j_Time];
  sprintf(Plot_Time, "Time = %5.2f", Current_Time);
  sprintf(Time_Eraser, "Time = %5.2f", Last_Time);
  cpgqch(&char_Size);
  cpgsch(2.0);
  cpgsci(0);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 0.0, Time_Eraser);
  cpgsci(1);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 0.5, Plot_Time);
  cpgsch(char_Size);
  free(Plot_Time); free(Time_Eraser);
  /*     E N D : Plotting Current Time on Top */

  if(Table->TYPE_of_INITIAL_CONDITION == 0 && Table->TYPE_of_MODEL == 0)
    assert( n == Table->No_of_INDIVIDUALS );

  free(x_Data);
  free(y_Data);
}

void Community_Scatter_Plot_Representation( Parameter_Table * Table,
					                                  int i_Replicate, int j_Time )
{
  int i, j, k;
  int N, n;
  double x, y;
  Parameter_CPGPLOT * CPG;
  int Sp;
  float Range_x[2], Range_y[2];
  int color_Index;
  int type_of_Line;
  int type_of_Width;
  int type_of_Symbol;

  Community * Patch; 
  Community ** P = Table->Patch_System;
  double * Y = Table->Vector_Model_Variables;

  N =  Total_Population(Y, Table);

  if (Table->TYPE_of_MODEL == 0 )
    assert( N == (Table->No_of_RESOURCES * Table->No_of_INDIVIDUALS) );

  if( N > 0 ) { /* Otherwise nothing to represent!!! */
    
    float * x_Data = (float *)calloc( N, sizeof(float) );
    float * y_Data = (float *)calloc( N, sizeof(float) );

    char ** Y_label = (char **)malloc( sizeof(char *) * Table->SUB_OUTPUT_VARIABLES );
    for(i=0; i < Table->SUB_OUTPUT_VARIABLES; i++){
      k = Table->OUTPUT_VARIABLE_INDEX[i];
      // Y_label[i]  = P->Output_Variable_Name[k];
      Y_label[i]  = Table->Output_Variable_Symbol[k];
    }

    Range_x[0] = 0.0;  Range_x[1] = P[0]->X_DIMENSION;
    Range_y[0] = 0.0;  Range_y[1] = P[0]->Y_DIMENSION;

    type_of_Line = 1;
    type_of_Width = 1;
    cpgsls(type_of_Line);
    cpgslw(type_of_Width);
    cpgask( 0 );
    cpgslct(Table->CPG_STO->DEVICE_NUMBER);      /* Selecting Device */
    
    for(Sp = 0; Sp<Table->LOCAL_STATE_VARIABLES; Sp++) {
      
      n=0;
      for(i = 0; i<Table->No_of_CELLS; i++) {
	      Patch = P[i]; 
	      for(j = 0; j<Patch->n[Sp]; j++) {
	        x_Data[n] = (float)(Patch->center.x - 0.5 + gsl_rng_uniform(r));
	        /* This is because STEP is 1.0 */
	        y_Data[n] = (float)(Patch->center.y - 0.5 + gsl_rng_uniform(r));
	        /* This is because STEP is 1.0 */
	        n++;
	      }
      }
      
      if (Sp == 0) {
        color_Index    = 2;          /* 2: Red        */ 
	      type_of_Symbol = 1;          /* 1: .  (point) */
	      
        cpg_XY_scattered(n, x_Data, y_Data, Range_x, Range_y,
		                     color_Index, type_of_Symbol,
		                     "X", "Y", Y_label[Sp]);
      }
      else if (Sp == 1) {
        color_Index    = 7;          /* 7: Yellow     */
	      type_of_Symbol = 2;          /* 2:  +         */
	      
        cpg_XY_scattered(n, x_Data, y_Data,  Range_x, Range_y,
			                   color_Index, type_of_Symbol, 
                         "X", "Y", Y_label[Sp]);
      }
      else if (Sp == 2){
        color_Index    = 5;    /* 4: Blue  5: Light Blue  6: Magenta    */
	      type_of_Symbol = 3;    /* 2: +     3: *           4: o (cercle) */

	      cpg_XY_scattered(n, x_Data, y_Data,  Range_x, Range_y,
			                   color_Index, type_of_Symbol, 
                         "X", "Y", Y_label[Sp]);
      }
      else if (Sp == 3){
	      color_Index    = 6;    /* 4: Blue  5: Light Blue  6: Magenta    */
	      type_of_Symbol = 4;    /* 2: +     3: *           4: o (cercle) */ 
	      
        cpg_XY_scattered(n, x_Data, y_Data,  Range_x, Range_y,
			                   color_Index, type_of_Symbol, 
                         "X", "Y", Y_label[Sp]);
      }
      else assert(Table->SUB_OUTPUT_VARIABLES == 4);  

      cpgsls(1);
      cpgslw(1);
      /* n: Total number of individuals across species */
      // Print_Meta_Community_Patch_System (Table);
    
      /* B E G I N : Plotting Current Time on Top */
      char * Plot_Time  = (char *)calloc( 50, sizeof(char));
      char * Time_Eraser = (char *)calloc(50, sizeof(char));
      float char_Size;
      float x_Time_Position = 0.02 * Range_x[1]; 
      float y_Time_Position = 1.09 * Range_y[1];
      static double Current_Time  = 0.0; 
      double Last_Time            = Current_Time;
      Current_Time  = Table->T->Time_Vector[j_Time];
      sprintf(Plot_Time, "Time = %5.2f", Current_Time);
      sprintf(Time_Eraser, "Time = %5.2f", Last_Time);
      cpgqch(&char_Size);
      cpgsch(2.0);
      cpgsci(0);
      cpgptxt(x_Time_Position, y_Time_Position, 0.0, 0.5, Time_Eraser);
      cpgsci(1);
      cpgptxt(x_Time_Position, y_Time_Position, 0.0, 0.5, Plot_Time);
      cpgsch(char_Size);
      free(Plot_Time); free(Time_Eraser);
      /*     E N D : Plotting Current Time on Top */
    }

    if(Table->TYPE_of_INITIAL_CONDITION == 0 && Table->TYPE_of_MODEL == 0)
      assert( n == Table->No_of_INDIVIDUALS );

    free(Y_label);  
    free(x_Data);
    free(y_Data);
  }
}
