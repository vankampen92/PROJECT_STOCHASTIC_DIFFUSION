#include <MODEL.h>
#include "IO_Procedures_AIDS.h"

void Reading_Demographic_Parameters_from_File(char * File_Name, double ** Data, int * N,
					      int No_of_PARAMETERS) 
{
  int j,n; 
  char * Dummy = (char *)calloc( 50, sizeof(char) );
  double y; 
  FILE * f;

  printf("\n [Entering function Reading_Matrix_from_File(...)]\n");
  printf(" Reading File %s...\n", File_Name);
  if((f=fopen(File_Name,"r")) == NULL) {
    printf("File non-existent! Cannot open file.\n");
    printf("Program aborted!!!"); exit(1);
  }
  

  for(j=0; j<No_of_PARAMETERS + 3; j++)
    if(j ==  (No_of_PARAMETERS+2)  ) fscanf(f, "%s\n", Dummy);
    else                             fscanf(f, "%s\t", Dummy);

  n=0;
  while ( fscanf(f, "%s\t", Dummy) != EOF ){
  
    for(j=0; j < No_of_PARAMETERS; j++) {

      fscanf(f, "%lf\t", &y);
      Data[n][j] = y; 
    }
    
    fscanf(f, "%lf\t%lf\n", &y, &y);
    
    n++;
  }
 
  * N = n;  //Number of Different Parameter sets  

  fclose(f);
  free(Dummy);
  
  printf(" File %s has been read successfully\n", File_Name);
  printf(" [Exiting function Reading_Matrix_from_File(...)]\n\n");
}

void Writing_Demographic_Parameters_Matrix(Parameter_Table * Table, const char * City_Name,
					   double ** Data,
					   int N,
					   int * Index, int No_of_PARAMETERS) 
{
  int j,n; 
  int key;

  printf("%s\t", "City");

  for(j=0; j < No_of_PARAMETERS; j++) {
    key = Index[j];
    printf("%s\t", Table->Symbol_Parameters[key]);
  }
  printf("\n");
      
  n=0;
  while ( n < N ){

    printf("%s\t", City_Name);
    
    for(j=0; j < No_of_PARAMETERS; j++) printf("%g\t", Data[n][j]);

    printf("\n"); 
    
    n++;
  }
}
