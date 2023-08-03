#include <MODEL.h>

void Observed_Data_Alloc( Observed_Data * Data, int No_of_OBSERVED_VARIABLES, int No_of_TIMES)
{                                                   
  int i;
  
  Data->N = (double **)calloc( No_of_OBSERVED_VARIABLES, sizeof(double *) );
  for(i = 0; i<No_of_OBSERVED_VARIABLES; i++)
    Data->N[i] = (double *)calloc( No_of_TIMES, sizeof(double) );

  Data->Name = (char *)calloc(100, sizeof(char)); 

  Data->Time_Data_Vector = (double *)calloc( No_of_TIMES, sizeof(double) );
}

void Observed_Data_Initialization( Observed_Data * Data, int No_of_OBSERVED_VARIABLES,
				                           int No_of_TIMES, double ** Data_Matrix,
				                           const char * Name )
{
  int i, j;
  
  Data->No_of_VARIABLES = No_of_OBSERVED_VARIABLES;
  Data->No_of_POINTS    = No_of_TIMES;

  for(i = 0; i<No_of_OBSERVED_VARIABLES; i++)
    for(j = 0; j<No_of_TIMES; j++)
      Data->N[i][j] = Data_Matrix[i][j];

  strcpy(Data->Name, Name); 
}

void Observed_Data_Free( Observed_Data * Data)
{
  int i;
  
  for(i = 0; i<Data->No_of_VARIABLES; i++) free(Data->N[i]);
  free(Data->N);
  
  free(Data->Name);

  free(Data->Time_Data_Vector);  
}
