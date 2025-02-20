#include <MODEL.h>

void Printing_Putative_Recipient_Lists(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 
  int i_List, k_List; 

  for(i = 0; i<Table->No_of_RESOURCES; i++) {
   Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 

    printf("\n\n");
    printf("Putative Recipient List of Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", i, i_Focal, k_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]:\n");   

    for(j=0; j<Table->Putative_Recipient_List_Indeces[i][Table->No_of_PROFILES]; j++) {

      Calculate_Strain_and_Profile(Table, Table->Putative_Recipient_List_Indeces[i][j], &i_List, &k_List);   

      printf("Strain ID %d:\t Strain TYPE (%d) and Profile (%d): [ ", 
        Table->Putative_Recipient_List_Indeces[i][j], i_List, k_List);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_List][k_List][k]);
      printf("]\n"); 
    }

    printf("\n");
    printf("Putatative Recipient List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Putative_Recipient_List_Indeces[i][Table->No_of_PROFILES]; j++)
      printf(" %d ", Table->Putative_Recipient_List_Indeces[i][j]);
    printf("}");

    printf("\n\n");
  }
}
