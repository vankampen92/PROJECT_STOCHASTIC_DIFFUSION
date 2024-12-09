#include <MODEL.h>

void Print_Strain_Profile(int * Profile, int No_of_PLASMIDS) 
{
    /* This input Profile should have been created to store
       (No_of_PLASMIDS + 1) integers 
    */
      int k;
      printf("[ ");
      for(k=0; k<No_of_PLASMIDS; k++) 
        printf("%d ", Profile[k]);
      printf("] ");
      printf("[[ %d ]]", Profile[No_of_PLASMIDS]);   
}