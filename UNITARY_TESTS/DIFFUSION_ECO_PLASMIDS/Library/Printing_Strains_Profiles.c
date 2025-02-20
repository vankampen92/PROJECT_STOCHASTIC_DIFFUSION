#include <MODEL.h>

void Printing_Strains_Profiles(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 

  for(i=0; i<Table->No_of_STRAINS; i++)
    printf("Strain Type: %d\t No of Profiles of this strain type: %d\n", i, Table->n[i]);
  printf("\n"); printf("\n");

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
    Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 
   
    printf("Strain ID: %d\t Strain Type: %d\t Strain Profile [ ", i, i_Focal);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]");

    printf("\n");  
  }
}
