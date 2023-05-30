#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* Compilation:
 * :~$ gcc -c main.c 
 * 
 * Link: 
 * :~$ gcc -o test main.o -lm
 *
 * Execution
 * :~$ ./test
 */
 
unsigned long long int Fun_Hyper(unsigned int D,
				 unsigned int k)
{
  /* Number of configuations... */
  unsigned long long int S, i; 
  
  if(D == 0) return 1; 

  S = 0; 
  for(i=0; i<=k; i++)
    S += Fun_Hyper(D-1, i);  
  
  return S;
}

void Groups(unsigned int D, int D_Max, int * V, int * s)
{
  /* Vector V[] should be declared with enough room to host all values */
  int i, j, k; 
  int n_0, n; 

  if (D < D_Max ) {
    n_0 = V[0]; 

    j=0;
    k=0; 
    while ( n_0 > 0 )
    {
      for(i=0; i < *s; i++) {
        n = V[i]-k;
        if(n > 0)  
          V[j++] = n;  
      }
      n_0--;
      k++; 
    }
    *s = j;  
    Groups(D+1, D_Max, V, s);
  }
}

int  main()
{  
   int * V; 
   int i, j, s_Max; 
   unsigned long long int S, SBF;
   int s;

   unsigned int D   = 2; 
   unsigned int n_0 = 5;
   
   printf("Enter the Number of Resources (Dimensions)... ");
   scanf("%d", &D); 
   printf("Number of Resources (Dimensions): %d\n", D);
   
   printf("Enter the Number of Consumers (Total)... ");
   scanf("%d", &n_0); 
   printf("Number of Consumers (Total): %d\n", n_0);
   
   getchar();
   
   S = Fun_Hyper(D, n_0);
   printf("Number of configurations S(D, n_0) = S(%d, %d) = %lld\n\n", D, n_0, S);

   printf("Vector of Sub-groups per Dimension: \n");
   for( j=1; j<=D; j++) {

    s_Max = Fun_Hyper(j-1, n_0);
    V = (int *)calloc(s_Max, sizeof(int));

    V[0] = n_0+1; 
    s = 1; 
    Groups(1, j, V, &s);

    printf("V(D = %d) = [ ", j);
    for(i=0; i < s; i++) 
      printf(" %d ", V[i]);
    
    printf(" ]\n");

    free(V);
   }

   return 0;
}
