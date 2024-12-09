#include <MODEL.h>

int SumUP_Profile ( int * Profile, int N )
{
  int i, S;

  S = 0; 
  for( i=0; i<N; i++)
    S += Profile[i];

  return(S);  
}