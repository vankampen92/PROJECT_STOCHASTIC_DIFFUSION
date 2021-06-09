#include <MODEL.h>
#include "IO_Procedures_Aux.h"

void Reading_Observed_Data(char * OBSERVED_DATA_FILE, double ** Data_Matrix,
			   int No_of_OBSERVED_VARIABLES_MAX, int No_of_COLUMNS,
			   int * No_of_OBSERVED_VARIABLES,
			   int * Output_Variable_Index, int N,
			   double * Initial_Observed_Adult_Population,
			   double * Initial_Observed_Female_Sexual_Workers,
			   double * Initial_Observed_Female_Sexual_Workers_Infected )
{
  int k, i, n, m;
  int key;
  
  char ** Name_of_Rows;
  double * Time_Vector;
  int Index[3] = {6, 4, 14};  /* These are the potential set of output variables that may 
				 have been also observed */
  
  int * Last_Output_Variable_Index       = (int *)calloc( N-3, sizeof(int) ); 
  int * Index_Output_Variable_Observed   = (int *)calloc( 3, sizeof(int) );
  
  double ** Data = (double **)calloc(No_of_OBSERVED_VARIABLES_MAX, sizeof(double *) ) ;
  for (i=0; i<No_of_OBSERVED_VARIABLES_MAX; i++)
    Data[i] = (double *)calloc(No_of_COLUMNS, sizeof(double) );
    
  Reading_Standard_Data_Matrix_from_File( OBSERVED_DATA_FILE,
					  Data,
					  &No_of_OBSERVED_VARIABLES_MAX, No_of_COLUMNS, 
					  0, Name_of_Rows,
					  0, Time_Vector );

  * Initial_Observed_Adult_Population               = Data[0][0];
  * Initial_Observed_Female_Sexual_Workers          = Data[1][0];
  * Initial_Observed_Female_Sexual_Workers_Infected = Data[2][0]; 
  
  assert( No_of_OBSERVED_VARIABLES_MAX == 3);

  bool Bool = isDisjoint(Index, 3, Output_Variable_Index, N);

  if (Bool) (* No_of_OBSERVED_VARIABLES) = 0;
  else {
    Bool = 1; 
    
    n=0; 
    for( i=3; i<N; i++ ) Last_Output_Variable_Index[n++] = Output_Variable_Index[i];

    if(N > 3) Bool = isDisjoint(Index, 3, Last_Output_Variable_Index, N-3);

    if(Bool) { 
      m = 0; n=1; 
      while( m<3 || n==1 ) {
	n = 0; 
	if( Output_Variable_Index[m] == Index[0] ) { key = 0; n = 1; }
	if( Output_Variable_Index[m] == Index[1] ) { key = 1; n = 1; } 
	if( Output_Variable_Index[m] == Index[2] ) { key = 2; n = 1; } 
	
	if( n == 1 ) {
	  
	  for (k=0; k<No_of_COLUMNS; k++)  Data_Matrix[m][k] = Data[key][k];
	  
	  m++; 
	}	
	
      }

      * No_of_OBSERVED_VARIABLES = m; 
      
      if( m<3 ) { 
	for(k=m; k<3; k++) Index_Output_Variable_Observed[k-m] = Output_Variable_Index[k];
      
	Bool = isDisjoint(Index, 3, Index_Output_Variable_Observed, 3-m);

	if( !Bool ){
	  printf(" Wrong Output Variable List!!!\n");
	  printf(" Output variables 6, 4 and 14 do not appear in the right order\n"); 
	  printf(" in the three first positions of the Output_Variable_Index vector\n");
	  printf(" At least of the these three frist positions is occupied by a non-observed\n");
	  printf(" output variable\n");
	  printf(" The whole Output_Variable_Index vector is:\n");
	  printf(" [ "); 
	  for( i=0; i<N; i++ ) printf(" %d ", Output_Variable_Index[i]);
	  printf(" ] ");
	  printf(" The program will exit\n");
	  exit(0);
	}
      }
    }
    else {
      printf(" Wrong Output Variable List!!!\n");
      printf(" Output variables 6, 4 and 14 appear after the third position in\n");
      printf(" the Output_Variable_Index vector (at least one of them)\n");
      printf(" The whole Output_Variable_Index vector is:\n");
      printf(" [ "); 
      for( i=0; i<N; i++ ) printf(" %d ", Output_Variable_Index[i]);
      printf(" ] ");
      printf(" The program will exit\n");
      exit(0);
    }
  }

  free(Last_Output_Variable_Index);    
  free(Index_Output_Variable_Observed);

  for (i=0; i<No_of_OBSERVED_VARIABLES_MAX; i++) free(Data[i]); 
  free(Data);
}
 
