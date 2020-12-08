#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Community_Scatter_Plot_Representation( Parameter_Table * Table, int SAME_PLOT, int Sp )
{
  int i,j;
  int N, n;
  double x, y; 
  
  Community ** P = Table->Patch_System;
  double * Y = Table->Vector_Model_Variables;
  
  N =  Total_Population(Y, Table);

  if (Table->TYPE_of_MODEL == 0 )
    assert( N == Table->No_of_INDIVIDUALS );

  double * x_Data = (double *)calloc( N, sizeof(double) );
  double * y_Data = (double *)calloc( N, sizeof(double) );

  n=0; 
  for(i = 0; i<Table->No_of_CELLS; i++) 
    for(j = 0; j<P[i]->n[Sp]; j++) {
      x_Data[n] = P[i]->center.x - 0.5 + gsl_rng_uniform(r); /* This is because STEP is 1.0 */
      y_Data[n] = P[i]->center.y - 0.5 + gsl_rng_uniform(r); /* This is because STEP is 1.0 */
      n++; 
    }
  
  assert( n == N ); 
  
  Table->CPG->CPG_RANGE_X_0 = 0.0;  Table->CPG->CPG_RANGE_X_1 = P[0]->X_DIMENSION;
  Table->CPG->CPG_RANGE_Y_0 = 0.0;  Table->CPG->CPG_RANGE_Y_1 = P[0]->Y_DIMENSION;

  Table->CPG->color_Index  = 1 + Sp;
  Table->CPG->type_of_Line = 2;
  Table->CPG->type_of_Width = 1;
  Table->CPG->type_of_Symbol = 1; 

  CPGPLOT___X_Y___S_C_A_T_T_E_R_E_D___S_A_M_E___P_L_O_T ( Table->CPG,
							  SAME_PLOT,
							  N, 
							  x_Data, y_Data,
							  "", 
							  "", 
							  "",
							  1, 1 ); 
  free(x_Data);
  free(y_Data); 
}
