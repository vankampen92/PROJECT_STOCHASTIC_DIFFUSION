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

void Plotting_Current_Time_On_Top(float *Range_x, float *Range_y, Parameter_Table * Table, int j_Time)
{
  char *Plot_Time = (char *)calloc(50, sizeof(char));
  char *Time_Eraser = (char *)calloc(50, sizeof(char));
  float char_Size;
  float x_Time_Position = 0.5 * Range_x[1];
  float y_Time_Position = 1.09 * Range_y[1];
  static double Current_Time = 0.0;
  double Last_Time = Current_Time;
  Current_Time = Table->T->Time_Vector[j_Time];
  sprintf(Time_Eraser, "Time = %5.2f", Last_Time);
  sprintf(Plot_Time, "Time = %5.2f", Current_Time);
  cpgqch(&char_Size);
  cpgsch(2.0);
  cpgsci(0);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 0.5, Time_Eraser);
  cpgsci(1);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 0.5, Plot_Time);
  cpgsch(char_Size);
  free(Plot_Time);
  free(Time_Eraser);
}

void Community_Bar_Plot_Representation(Parameter_Table * Table, int j_Time) 
{  
  int SAME_PLOT = 0; 
  int i,j,k;
  int R, RP; 

  assert( Table->TYPE_of_MODEL == 20); /* MODEL = DIFFUSION_ECOEVO_PLANTS */

  double * X = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));
  double * Y = (double *)calloc(Table->No_of_RESOURCES, sizeof(double)); /* Every entry to zero!!! */
 
  for(i=0; i<Table->No_of_RESOURCES; i++)  X[i] = (double)i+1.0; /* Species labeled 
                                                                    as 1 to No_of_RESOURCES types 
                                                                  */
  for(k = 0; k<Table->No_of_RESOURCES; k++) {  
    for (j=0; j<Table->No_of_CELLS; j++) {
    
      RP  = j*Table->LOCAL_STATE_VARIABLES + 2*k;
      R   = j*Table->LOCAL_STATE_VARIABLES + 2*k+1;

      Y[k] += Table->Vector_Model_Variables[R];  /* Established adult plants across cells */
    }
  } 
  
  if (j_Time == 0) SAME_PLOT = 0; 
  else             SAME_PLOT = 1;  /* Always redrawing the plot.
                                      Otherwise, SAME_PLOT = 1, for j_Time != 0 */ 

  float * Range_x = (float *)calloc(2, sizeof(float));
  float * Range_y = (float *)calloc(2, sizeof(float));

  int SCALE_X = 1; /* Respect always this given rage for the X axis  */
  Table->CPG->CPG_RANGE_X_0 = 0.0;   
  Table->CPG->CPG_RANGE_X_1 = (double)Table->No_of_RESOURCES + 1.0;
  /* Y axis range is controlled by command line 1 (according to carrying capacity, i.e., -HK 1000) */
  Range_x[0] = Table->CPG->CPG_RANGE_X_0; Range_x[1] = Table->CPG->CPG_RANGE_X_1;
  Range_y[0] = Table->CPG->CPG_RANGE_Y_0; Range_y[1] = Table->CPG->CPG_RANGE_Y_1; 
                     
  /* B E G I N : Plotting Current Time on Top */
  // Plotting_Current_Time_On_Top(Range_x, Range_y, Table, j_Time);
  /*     E N D : Plotting Current Time on Top */

  /* Deleting the title: Table->CPG->Title[0] = '\0'; */
  Table->CPG->color_Index = 1+j_Time; 
  CPGPLOT___B_A_R___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( Table->CPG,
                                                          SAME_PLOT,
                                                          Table->No_of_RESOURCES,
                                                          X, Y,
                                                          "Species Phenotypes",
                                                          "Abundance",
                                                          "", 
                                                          SCALE_X, Table->CPG->CPG_SCALE_Y );

  free(Range_x); free(Range_y); 
  free(X); free(Y);
}