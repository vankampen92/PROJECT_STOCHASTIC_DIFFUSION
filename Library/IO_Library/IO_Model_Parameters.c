#include <MODEL.h>
#include "IO_Procedures_Aux.h"

void Reading_Model_Parameters_from_File(char * File_Name, double ** Data, int * N,
					int No_of_PARAMETERS) 
{
  int j,n; 
  char * Dummy = (char *)calloc( 50, sizeof(char) );
  double y; 
  FILE * f;

  printf("\n [Entering function Reading_Model_Parameters_from_File(...)]\n");
  printf(" Reading File %s...\n", File_Name);
  if((f=fopen(File_Name,"r")) == NULL) {
    printf("File non-existent! Cannot open file.\n");
    printf("Program aborted!!!"); exit(1);
  }

  /* Reading the Title Columns */
  for(j=0; j<No_of_PARAMETERS + 1; j++)
    if(j == (No_of_PARAMETERS) ) fscanf(f, "%s\n", Dummy);
    else                         fscanf(f, "%s\t", Dummy);

  n=0;
  while ( fscanf(f, "%lf\t", &y) != EOF ){
 
    Data[n][0] = y;   
    
    for(j=1; j < No_of_PARAMETERS; j++) {

      fscanf(f, "%lf\t", &y);
      Data[n][j] = y; 
    }
    
    fscanf(f, "%lf\n", &y);
    Data[n][j] = y; 
    
    n++;
  }
 
  * N = n;  //Number of Different Parameter sets  

  fclose(f);
  free(Dummy);
  
  printf(" File %s has been read successfully\n", File_Name);
  printf(" [Exiting function Reading_Model_Parameters_into_File(...)]\n\n");
}

void Writing_Model_Parameters_into_File(char * File_Name, char ** Title_Row, double ** Data,
					int N, int No_of_PARAMETERS) 
{
  int i, j; 
  FILE * f;

  printf("\n [Entering function Writing_Model_Parameters_into_File(...)]\n");
  printf(" Writing File %s...\n", File_Name);
  f=fopen(File_Name,"w");

  for(j=0; j<No_of_PARAMETERS + 1; j++)
    if(j == (No_of_PARAMETERS) ) fprintf( f, "%s\n", Title_Row[j] );
    else                         fprintf( f, "%s\t", Title_Row[j]);
  
  for( i=0; i<N; i++ ) {   

    for(j=0; j<No_of_PARAMETERS+1; j++)
      if(j == (No_of_PARAMETERS) ) fprintf(f, "%g\n", Data[i][j]); 
      else                         fprintf(f, "%g\t", Data[i][j]); 
  }
  
  fclose(f);
  
  printf(" File %s has been written successfully\n", File_Name);
  printf(" [Exiting function Writing_Model_Parameters_into_File(...)]\n\n");
}

void Writing_Model_Parameters_Matrix(Parameter_Table * Table, const char * Label,
				     double ** Data,
				     int N, int No_of_PARAMETERS) 
{
  int i, j,n; 

  int TOTAL_No_of_Fitting_Parameters;
  TOTAL_No_of_Fitting_Parameters = Table->TOTAL_No_of_MODEL_PARAMETERS + Table->No_of_ERROR_PARAMETERS + Table->No_of_IC;

  assert(No_of_PARAMETERS == TOTAL_No_of_Fitting_Parameters);
  
  char ** Title_Parameters = (char **)calloc(TOTAL_No_of_Fitting_Parameters+1,
					     sizeof(char *) );
  for(i=0; i < TOTAL_No_of_Fitting_Parameters+1; i++)
    Title_Parameters[i] = (char *)calloc(100, sizeof(char));
  
  Creating_Title_Row (Table, Title_Parameters);
  
  n=0;
  while ( n < N ){

    printf("%s %d:\t", Label, n);

    for(j=0; j < No_of_PARAMETERS + 1; j++) {
      if(j == (No_of_PARAMETERS) ) printf("%g\n", Data[n][j]); 
      else                         printf("%g\t", Data[n][j]);
    }

    n++;
    
    Print_Press_Key(1,0,"."); 
  }

  for(i=0; i < TOTAL_No_of_Fitting_Parameters+1; i++) free(Title_Parameters[i]);
  free(Title_Parameters);
}

void Creating_Title_Row (Parameter_Table * Table, char ** Title_Parameters)
{
  /* Creating Title Row */
  /* B E G I N :   Defining first row of "Full_Parameter_Set.dat" file
   */
  int i, j, key;
  char * pTitle;
  char * pValue = (char *)calloc(20, sizeof(char) );

  int TOTAL_No_of_Fitting_Parameters;
  TOTAL_No_of_Fitting_Parameters = Table->TOTAL_No_of_MODEL_PARAMETERS + Table->No_of_ERROR_PARAMETERS + Table->No_of_IC;
  
  for(i=0; i<Table->TOTAL_No_of_MODEL_PARAMETERS; i++) {
    key = Table->Index[i];

    Title_Parameters[i][0]='\0'; 
    pTitle = strcat(Title_Parameters[i], Table->Symbol_Parameters[key]);
  }
  for(i=0; i<Table->No_of_IC; i++) {
    key = Table->IC_Space->Parameter_Index[i];
    
    j = Table->TOTAL_No_of_MODEL_PARAMETERS + i; 
    pTitle = strcat(Title_Parameters[j], Table->Model_Variable_Symbol[key]);
  }
  for(i=0; i<Table->No_of_ERROR_PARAMETERS; i++) {
    key = Table->E_Space->Parameter_Index[i];

    j = Table->TOTAL_No_of_MODEL_PARAMETERS + Table->No_of_IC + i;
    pTitle = strcat(Title_Parameters[j], "Err(");
    pValue[0]='\0'; sprintf(pValue, "%d", i);
    pTitle = strcat(Title_Parameters[j], pValue);
    pTitle = strcat(Title_Parameters[j], ")");
  }

  j = Table->TOTAL_No_of_MODEL_PARAMETERS + Table->No_of_IC + i;
  assert( j == TOTAL_No_of_Fitting_Parameters);
  
  j = Table->TOTAL_No_of_MODEL_PARAMETERS + Table->No_of_IC + Table->No_of_ERROR_PARAMETERS;
  pTitle = strcat(Title_Parameters[j], "NegLogLike");
  /*     E N D : --------------------------------------------------------
   */

  free(pValue);
}
