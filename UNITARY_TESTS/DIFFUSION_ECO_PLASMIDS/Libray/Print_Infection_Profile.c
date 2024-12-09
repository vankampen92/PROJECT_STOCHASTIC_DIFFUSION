#include <MODEL.h>

void Print_Infection_Profile(int * Profile, int i_List, Parameter_Table * Table)
{
  int k;
      printf("[ ");
      
      for(k=0; k<Table->No_of_PLASMIDS; k++) 
        if(Profile[k] == 1 ) {
          if( Table->IBP[i_List][k] == 1.0) 
            printf("Y ");
          else
            printf("N ");
        }
        else
          printf("%d ", Profile[k]);

      printf("] ");  
}
