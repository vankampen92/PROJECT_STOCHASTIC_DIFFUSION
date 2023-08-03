#include <HEADERS.h>
#include "IO_Procedures_Standard.h"

void Reading_Standard_Data_Matrix_from_File( char * File_Name,
					     double ** Matrix_Data,
					     int * No_of_Rows, int No_of_Columns,
					     int NAMES, char ** Name_of_Rows,
					     int TIMES, double * Time_Vector )
{
  int n;
  int j;
  double y;
  char * Name_Dummy = (char *)calloc( 50, sizeof(char) );

  FILE * fp;
  printf("\n [Entering function Reading_Matrix_from_File(...)]\n");
  printf(" Reading File %s...\n", File_Name);
  if((fp=fopen(File_Name,"r")) == NULL) {
    printf("File non-existent! Cannot open file.\n");
    printf("Program aborted!!!"); exit(1);
  }

  if(TIMES == 1) {

    if(NAMES == 1) { 
      fscanf(fp, "%s\t", Name_Dummy );
      for(j=0; j < No_of_Columns; j++) {
	if( j == (No_of_Columns - 1) ) fscanf(fp, "%lf\n", &y);
	else                         fscanf(fp, "%lf\t", &y);
	Time_Vector[j] = y;
      }
    }
    else {
      for(j=0; j < No_of_Columns; j++) {
	if( j == (No_of_Columns - 1) ) fscanf(fp, "%lf\n", &y);
	else                         fscanf(fp, "%lf\t", &y);
	Time_Vector[j] = y;
      }
    }
  }

  if(TIMES == 2) {

    if(NAMES == 1) { 
      fscanf(fp, "%s\t", Name_Dummy );
      for(j=0; j < No_of_Columns; j++) {
	if( j == (No_of_Columns - 1) ) fscanf(fp, "%s\n", Name_Dummy);
	else                         fscanf(fp, "%s\t", Name_Dummy);
	// Names_of_Columns[j] = Name_Dummy; 
      }
    }
    else {
      for(j=0; j < No_of_Columns; j++) {
	if( j == (No_of_Columns - 1) ) fscanf(fp, "%s\n", Name_Dummy);
	else                         fscanf(fp, "%s\t", Name_Dummy);
	// Names_of_Columns[j] = Name_Dummy; 
      }
    }
  }
  
  
  if (NAMES == 1) {
    n=0;
    while ( fscanf(fp, "%s\t", Name_Dummy) != EOF ){
      Name_of_Rows[n][0] = '\0'; 
      char * p = strcat( Name_of_Rows[n], Name_Dummy ); 

      // strcpy(Name_of_Rows[n], Name_Dummy);
      
      for(j=0; j < No_of_Columns; j++) {
	if( j == (No_of_Columns - 1) ) fscanf(fp, "%lf\n", &y);
	else                           fscanf(fp, "%lf\t", &y);
	Matrix_Data[n][j] = y;
      }
      
      n++;
    }
    
    fclose(fp);
    (*No_of_Rows) = n;
  }
  else {
    n=0; 
    while ( fscanf(fp, "%lf\t", &y) != EOF ){
      
      Matrix_Data[n][0] = y;
      
      for(j=1; j < No_of_Columns; j++) {
	if( j == (No_of_Columns - 1) ) fscanf(fp, "%lf\n", &y);
	else                           fscanf(fp, "%lf\t", &y);
	Matrix_Data[n][j] = y;
      }
      
      n++;
    }
    
    fclose(fp);
    (*No_of_Rows) = n;
  }
  
  printf(" File %s has been read successfully\n", File_Name);
  printf(" [Exiting function Reading_Matrix_from_File(...)]\n\n");

  free(Name_Dummy);
  // printf("Press any key to continue...\n"); getchar();
}

void Writing_Standard_Data_Matrix(double ** Matrix_Data,
				  int N_row, int N_column,
				  int NAMES, char ** Name_of_Rows,
				  int TIMES, double * Time_Vector)

{
  /* Input parameters:
     . N_row: Number of No_of_Rows
     . N_column: Number of Sampling Times
     . Matrix_Data[][]: Data
     . Name_of_Rows[][]: No_of_Rows Names
     . Time_Vector[]: a double array containing sampling times
  */
  int i,j;

	if(TIMES == 1) {

	  if(NAMES == 1) fprintf(stdout, "%s\t", "Name of Rows");
	  
	  for(j=0; j < N_column; j++) fprintf(stdout, "%g\t", Time_Vector[j]);
	  
	  fprintf(stdout, "\n");
	}

	for(i=0; i < N_row; i++){
	  
	  if(NAMES == 1) fprintf(stdout, "%s\t", Name_of_Rows[i] );
	  
	  for(j=0; j < N_column; j++) fprintf(stdout, "%g  ", Matrix_Data[i][j]);
	  
	  fprintf(stdout, "\n");
	}
}

void Writing_Standard_Data_Matrix_to_File(char * File_Name,
					  double ** Matrix_Data,
					  int N_row, int N_column, 
					  int NAMES, char ** Name_of_Rows,
					  int TIMES, double * Time_Vector)

{
  /* Input parameters:
     . N_row: Number of No_of_Rows
     . N_column: Number of Sampling Times
     . Matrix_Data[][]: Data
     . Name_of_Rows[][]: No_of_Rows Names
     . Time_Vector[]: a double array containing sampling times
  */
  int i,j;

  FILE * FP = fopen(File_Name, "w");
  
  if(TIMES == 1) {

    if(NAMES == 1) fprintf(FP, "%s\t", "Name of Rows");

    for(j=0; j < N_column; j++) fprintf(FP, "%g\t", Time_Vector[j]);

    fprintf(FP, "\n");
  }
  
  for(i=0; i < N_row; i++){
    
    if(NAMES == 1) fprintf(FP, "%s\t", Name_of_Rows[i] );
    
    for(j=0; j < N_column; j++) fprintf(FP, "%g\t", Matrix_Data[i][j]);
    
    fprintf(FP, "\n");
  }
  
  fclose(FP);
}
