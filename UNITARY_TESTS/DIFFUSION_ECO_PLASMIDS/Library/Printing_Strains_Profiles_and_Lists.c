#include <MODEL.h>

void Printing_Strains_Profiles_and_Lists(Parameter_Table * Table)
{
  int i, j, k; 
  int i_Focal, k_Focal; 
  int i_List, k_List; 

  for(i=0; i<Table->No_of_STRAINS; i++)
    printf("Strain Type: %d\t No of Profiles of this strain type: %d\n", i, Table->n[i]);
  printf("\n");

  /* Adjacent Lists: */
  for(i = 0; i<Table->No_of_RESOURCES; i++) {
   Calculate_Strain_and_Profile(Table, i, &i_Focal, &k_Focal); 
   
   printf("Strain ID: %d\t Strain Type: %d\t Strain Profile [ ", i, i_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]");
  
    printf("\n");
    printf("Competition List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Competition_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Competition_List_Indeces[i][j]);
    printf("}"); 

    printf("\n");
    printf("Conjugation List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Conjugation_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Conjugation_List_Indeces[i][j]);
    printf("}");

    printf("\n\n");
    printf("Recipient List of Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", i, i_Focal, k_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]:\n");   
    for(j=0; j<Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES]; j++) {

      Calculate_Strain_and_Profile(Table, Table->Recipient_List_Indeces[i][j], &i_List, &k_List);      
      printf("Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", Table->Recipient_List_Indeces[i][j], i_List, k_List);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_List][k_List][k]);
      printf("]\n"); 
    }
    printf("\n");
    printf("Recipient List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Recipient_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Recipient_List_Indeces[i][j]);
    printf("}");

    printf("\n\n");

    printf("Donor List of Strain ID %d:\t Bacterial Type (%d) and Profile (%d): [ ", i, i_Focal, k_Focal);
    for(k=0; k<Table->No_of_PLASMIDS; k++) 
      printf("%d ", Table->Strain_Profiles[i_Focal][k_Focal][k]);
    printf("]:\n");  
    for(j=0; j<Table->Donor_List_Indeces[i][Table->No_of_RESOURCES]; j++) {
      
      Calculate_Strain_and_Profile(Table, Table->Donor_List_Indeces[i][j], &i_List, &k_List);      
      printf("Strain ID %d:\t Bacterial Type (%d) and Profle (%d): [ ", Table->Donor_List_Indeces[i][j], i_List, k_List);
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        printf("%d ", Table->Strain_Profiles[i_List][k_List][k]);
      printf("]\n");      
    }
    printf("\n");
    printf("Donor List [%d: (%d, %d)] = {", i, i_Focal, k_Focal);
    for(j=0; j<Table->Donor_List_Indeces[i][Table->No_of_RESOURCES]; j++)
      printf(" %d ", Table->Donor_List_Indeces[i][j]);
    printf("}");
    printf("\n\n");
  }
  
  printf("\n");
  printf("Plasmid Compatibility List:\n");
  for(i=0; i<Table->No_of_PLASMIDS; i++) {
    printf("Plasmid ID: %d\t Compatibility list: [ ", i);
    for(j=0; j<Table->Plasmid_Compatibility_Indeces[i][Table->No_of_PLASMIDS]; j++)
      printf("%d ", Table->Plasmid_Compatibility_Indeces[i][j]);
    printf("]\n");
  }        
}