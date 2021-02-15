#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Community_Scatter_Plot_Representation( Parameter_Table * Table,
					    int i_Replicate, int j_Time )
{
  int i,j;
  int N, n;
  double x, x_0, c_x, y, y_0, c_y;
  static int DEVICE_NUMBER;
  Parameter_CPGPLOT * CPG;
  int Sp;
  float Range_x[2], Range_y[2];
  int color_Index  = 2 + Sp;
  int type_of_Line = 1;
  int type_of_Width = 2;
  int type_of_Symbol = 1;

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

  if(i_Replicate == 0 && j_Time == 1) {
    DEVICE_NUMBER = cpgopen( "/XSERVE" );
    cpgsubp(1, 1);
    cpgask( 0 );
  }

  for(Sp = 0; Sp<Table->LOCAL_STATE_VARIABLES; Sp++) {
    
    switch (Sp)
      {
      case 0: 
	x_0 = 0.0;
	y_0 = 0.0;
	break; 
      case 1:
	x_0 = Range_x[1]/2.0;
	y_0 = 0.0;
	break;
      case 2:
	x_0 = 0.0;
	y_0 = Range_y[1]/2.0; 
	break;
      case 3:
	x_0 = Range_x[1]/2.0;
	y_0 = Range_y[1]/2.0; 
	break;
      default:
	printf(" This plotting funciton does not work if the\n");
	printf(" number of species is greater than 4!!!\n");
	printf(" Here, the species index is %d,\n", Sp);
	printf(" which means that the number of species is\n");
	printf(" at least %d\n", Sp+1);
	printf(" The program will exit\n");
	exit(0); 
      }
	  
    n=0;
    for(i = 0; i<Table->No_of_CELLS; i++) {
      Patch = P[i];
      c_x = x_0 + 0.5 * Patch->center.x;
      c_y = y_0 + 0.5 * Patch->center.y;
      
      for(j = 0; j<Patch->n[Sp]; j++) {
	
	x_Data[n] = (float)(c_x - 0.25 + gsl_rng_uniform(r));
	/* This is because STEP is 1.0 */
	y_Data[n] = (float)(c_y - 0.25 + gsl_rng_uniform(r));
	/* This is because STEP is 1.0 */
	n++;
      }
    }
      
    type_of_Line = 1;
    type_of_Width = 2;
    cpgsls(type_of_Line);
    cpgslw(type_of_Width);
    cpgask( 0 );
    
    cpgslct(DEVICE_NUMBER);      /* Selecting Device */
    
    if (Sp == 0) {
      color_Index    = 3;
      type_of_Symbol = 1;
      cpg_XY_scattered(n, x_Data, y_Data, Range_x, Range_y,
		       color_Index, type_of_Symbol,
		       "X", "Y", "");
    }
    else {
      color_Index    = 3 + Sp;
      type_of_Symbol = 1 + Sp;
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
  Current_Time  = Table->T->Time_Vector[j_Time];
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
