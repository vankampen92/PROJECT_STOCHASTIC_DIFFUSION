#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Community_Scatter_Plot_Representation( Parameter_Table * Table,
					    int i_Replicate, int j_Time )
{
  int i,j;
  int N, n;
  double x, y;
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

  for(Sp = 0; Sp<Table->No_of_RESOURCES; Sp++) {

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
    type_of_Width = 2;
    cpgsls(type_of_Line);
    cpgslw(type_of_Width);
    cpgask( 0 );

    cpgslct(DEVICE_NUMBER);      /* Selecting Device */

    if (Sp == 0) {
      color_Index    = 2;
      type_of_Symbol = 17;
      cpg_XY_scattered(n, x_Data, y_Data, Range_x, Range_y,
		       color_Index, type_of_Symbol,
		       "X", "Y", "");
    }
    else {
      color_Index    = 2 + Sp;
      type_of_Symbol = 1 + Sp;
      cpg_XY_same_scattered(n, x_Data, y_Data,
			    color_Index, type_of_Symbol);
    }
    cpgsls(1);
    cpgslw(1);
  }

  /* n: Total number of individuals across species */
  // Print_Meta_Community_Patch_System (Table);

  if(Table->TYPE_of_INITIAL_CONDITION == 0 && Table->TYPE_of_MODEL == 0)
    assert( n == Table->No_of_INDIVIDUALS );


  free(x_Data);
  free(y_Data);
}
