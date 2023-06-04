#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct configurationinfo
{
  int * n;    /* Pointer to a configuration of the type (n_1, ..., n_D) of dimension D */
  int * anUp; /* Vector of indeces to neighgoring configurations of dimension D + 1    */
  int * anDw; 
              /* anUp[D] = a, where a is the no of UP neighboring configurations, a <= D   */
              /* anUp[D] = a, where a is the no of DOWN neighboring configurations, a <= D */
}configuration;

/* Compilation:
 * :~$ gcc -c main_Writing_Configuration_Table.c 
 * 
 * Link: 
 * :~$ gcc -o configuration main_Writing_Configuration_Table.o -lm
 *
 * Execution:
 * :~$ ./configuration
 */

unsigned long long int Func_Hyper(unsigned int D, unsigned int k)
{
  /* Exhaustive count of the number of configuations.
     
     For instance, for D=3, configurations, (n_x, n_y, n_z), should fulfill 
     the condition: 
                             n_x + n_y + n_z <= k,

     where n_x, n_y, and n_z can take values between 0 and k.
  */
  unsigned long long int S, i; 
  
  if(D == 0) return 1; 

  S = 0; 
  for(i=0; i<=k; i++)
    S += Func_Hyper(D-1, i);  
  
  return S;
}

void Groups(unsigned int D, int D_Max, int * V, int * s)
{
  /* This function calculates group configurations. 

     Vector V[] should be declared with enough room to host all groups 
  */
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

unsigned long long int General_Configuration_Table(unsigned int D, unsigned int n_0, 
                                                   int ** Configuration_Table)
{
  /* A recursion to build exhaustively all possible configurations such that, for instance, 
     for D=3, configurations, (n_x, n_y, n_z), fulfill the condition:

                             n_x + n_y + n_z <= n_0,

     where n_x, n_y, and n_z can take values between 0 and n_0

    Input: 
      . D,    dimension. 
      . n_0,  condition: n_x + n_y + n_z <= n_0, (for instance, for D=3), 
      . Configuration_Table[][] has been allocated before calling this function

    Output: 
      . Configuraiton_Table[][] is filled with all possible configurations. 
  */
  int i, j, k, l, s, s_0, d, g_D;
  unsigned long long int S;
  int * G_c; 

  int ***Configuration  = (int ***)calloc(D, sizeof(int **));

  int * S_D             = (int *)calloc(D, sizeof(int));
  for(i=0; i<D; i++) 
    S_D[i] = Func_Hyper(i+1, n_0); // S_1(n_0), ..., S_D(n_0) //
  
  for(i=0; i<D; i++){
    Configuration[i] = (int **)calloc(S_D[i], sizeof(int *)); 
    for(j=0; j < S_D[i]; j++)
      Configuration[i][j] = (int *)calloc(i+1, sizeof(int));
  }

  /* Index i=0 in Configuration[i][][] corresponds to the initial configuration, d=1.   */
  /* Index i=1 in Configuration[i][][] corresponds to the configuration, d=2, and so on */
  /* Configuration for d=1 is set up here: */
  for(i=0; i<=n_0; i++) 
    Configuration[0][i][0] = i;

  /* After each round of the loop the next configuration, d+1, is calculated 
     based on the previous one */
  d=1; 
  while (d < D)
  {
    g_D = Func_Hyper(d-1, n_0);
    G_c = (int *)calloc(g_D, sizeof(int));
    G_c[0] = n_0+1; 
    s = 1; 
    Groups(1, d, G_c, &s);

    S=0; 
    for(i=0; i<=n_0; i++) {
      s_0=0;
      for(k=0; k<s ; k++) {
        for(j=0; j<G_c[k]-i; j++) {
            for(l=0; l < d; l++){ 
              Configuration[d][S][l] = Configuration[d-1][j+s_0][l];
            }
            Configuration[d][S][d] = i; 
            S++;
        }
        s_0 += G_c[k];  
      }
    }  

    free(G_c);
    assert(S_D[d] == S);
    d++;
  }

  /* Output: Storing Configuration Table at D dimensions: */
  for(i=0; i < S; i++) 
     for(j=0; j < D; j++)
       Configuration_Table[i][j] = Configuration[D-1][i][j];

  for(i=0; i<D; i++){ 
    for(j=0; j < S_D[i]; j++)
      free(Configuration[i][j]);
    free(Configuration[i]);
  }
  free(Configuration);

  free(S_D);
  return(S);
}

/* Mapping of configuration states (n_1,..., n_D) into just one index (i) and viceversa */ 
int Configuration_to_i_Map(int * Configuration, int D, int n_0)
{
  int i, k, ks, n, m; 

  i = 0; ks=0;
  for (n=D-1; n >= 0; n--) {
    k  = Configuration[n];

    for (m = 0; m < k; m++)
      i += Func_Hyper(n, n_0-ks-m);
    
    ks+= k;
  }

  return i;
}

void i_to_Configuration_Map(int ** Configuration_Table, int i, int * Configuration, int D)
{
  int j;
  for(j=0; j < D; j++)
       Configuration[j] = Configuration_Table[i][j];
}

void Writing_General_Configuration_Table(int D,
                                         unsigned long long int S, 
                                         int ** Configuration_Table)
{
  /* S: Number of configuations... */
  int i, j;
  
  for(i=0; i < S; i++) {
    printf("Configuration[%d] = [", i);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table[i][j]);
    printf("]\n");
    if (i < 20 || i > S-20) getchar();
  }      
}      

void Comparing_Two_Configuration_Tables(int D, unsigned long long int S, 
                                        int ** Configuration_Table_D5, 
                                        int ** Configuration_Table)
{
  int i, j; 

  printf("Comparting two Configuration Tables...\n");
  
  for(i=0; i < S; i++) {
    printf("Configuration[%d] = [", i);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table_D5[i][j]);
    printf("]\t");

    printf("Configuration[%d] = [", i);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table[i][j]);
    printf("]\n");
    if (i < 20 || i > S-20) getchar();
  }
}

void Pickup_a_Configuration_from_Table(int D,
                                       int Index, 
                                       int ** Configuration_Table)
{
  /* Index: Order for a configuration... */
  int j;
  
  printf("Configuration[%d] = [", Index);
    for(j=0; j < D; j++)
      printf(" %d ", Configuration_Table[Index][j]);
  printf("]\n");
  // getchar();
}        

int Vector_of_Adjacent_Configurations(int * Configuration, int D, int n_0, int d,  
                                      int ** Nei, int * Index_Nei)
{
  int i, j, n, a; 

  n=0; 
  for(i=0; i<D; i++) 
    n += Configuration[i];   /* Controling Hyperplane Condtion!!! */
  
  a = 0;
  for(i=0; i<D; i++) {
     
    if( (Configuration[i]+d) >= 0 && (n+d) <= n_0 ) {
      for(j=0; j<D; j++) {
        if (j == i) Nei[a][j] = Configuration[j] + d;
        else        Nei[a][j] = Configuration[j];      
      }  
      a++;
    }  
  }

  for(i=0; i<a; i++) {
    Index_Nei[i] = Configuration_to_i_Map(Nei[i], D, n_0);
  }
  
  /* No of Neighboring configurations */
  return(a);
}

void Writing_Configuration_Table(unsigned int D, unsigned int n_0, int ** Configuration_Table);
unsigned long long int Func_Hyper_Brut_Force_D2(unsigned int D, unsigned int n_0);
unsigned long long int Func_Hyper_Brut_Force_D5(unsigned int D, unsigned int n_0);
void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0); 
int main()
{
  int i, k;
  unsigned long long int S, SBF;
  unsigned long long int No_of_CONFIGURATIONS;

  unsigned int D = 2;
  unsigned int n_0 = 5;

  printf("Enter the Number of Resources (Dimensions)... ");
  scanf("%d", &D);
  printf("Number of Resources (Dimensions): %d\n\n", D);
  // getchar();

  printf("Enter the Number of Consumers (Total)... ");
  scanf("%d", &n_0);
  printf("Number of Consumers (Total):%d\n\n", n_0);
  // getchar();

  S = Func_Hyper(D, n_0);
  printf("(1) Number of configurations S(D, n_0) = S(%d, %d) = %lld\n\n", D, n_0, S);

  int **Configuration_Table = (int **)calloc(S, sizeof(int *));
  for (i = 0; i < S; i++)
    Configuration_Table[i] = (int *)calloc(D, sizeof(int));
 
  int **Configuration_Table_D5 = (int **)calloc(S, sizeof(int *));
  for (i = 0; i < S; i++)
    Configuration_Table_D5[i] = (int *)calloc(5, sizeof(int));

  S = General_Configuration_Table(D, n_0, Configuration_Table);
  printf("(2) Number of configurations S(D, n_0) = S(%d, %d) = %lld\n\n", D, n_0, S);

  // Brut Force Calculation as a check
  if (D == 2 || D == 5)
  {
    if (D == 2)
      SBF = Func_Hyper_Brut_Force_D2(D, n_0);
    if (D == 5)
      SBF = Func_Hyper_Brut_Force_D5(D, n_0);

    printf("Brut force calculation (as a check) for dimensions 2 or 5...\n");
    printf("Number of configurations S(D, n_0) = S(%d, %d) = %lld\n\n", D, n_0, SBF);
    getchar();
  }

  if (D == 5) {
    Writing_Configuration_Table(D, n_0, Configuration_Table_D5);
    Comparing_Two_Configuration_Tables(D, S, Configuration_Table_D5, Configuration_Table);
  }
  printf("Writing General Configuration Table...\n");
  Writing_General_Configuration_Table(D, S, Configuration_Table);

  int ** Configuration = (int **)calloc(1, sizeof(int *));
  Configuration[0] = (int *)calloc(D, sizeof(int)); 
  printf("Enter a %d-dimensional configurational state... ", D);
  for (i=0; i<D; i++) {
    scanf("%d", &k);
    printf("Configuration[%d] = %d\n", i, k);
    Configuration[0][i] = k;  
  }
  printf("(1) You entered input configuration: \n");
  Writing_General_Configuration_Table(D, 1, Configuration);

  int Index = Configuration_to_i_Map(Configuration[0], D, n_0);
  printf("Configuration Index: %d\n", Index);

  i_to_Configuration_Map(Configuration_Table, Index, Configuration[0], D);
  printf("(2) Output Configuration: \t");
  Writing_General_Configuration_Table(D, 1, Configuration);
  Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
  
  /* There are up to D neighboring configurations around a given entered configuration */
  int a; /* No of Neighboring Configurations (a <= D)                                  */
  int * Index_Nei = (int *)calloc(D, sizeof(int));
  int ** Nei = (int **)calloc(D, sizeof(int *));
  for(i=0; i<D; i++)
    Nei[i] = (int *)calloc(D, sizeof(int));
  
  a = Vector_of_Adjacent_Configurations(Configuration[0], D, n_0, +1, Nei, Index_Nei);
  printf(" (1) There are %d neighboring configurations around ", a);
  Writing_General_Configuration_Table(D, 1, Configuration);  
  printf(" transitioning DOWN into it:\n");
  Writing_General_Configuration_Table(D, a, Nei);

  printf(" (2) There are %d neighforing configurations around ", a);
  Writing_General_Configuration_Table(D, 1, Configuration);
  printf(" transitioning DOWN into it:\n\n");  
  for(i=0; i<a; i++){
    Index = Index_Nei[i];
    Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
  }

  a = Vector_of_Adjacent_Configurations(Configuration[0], D, n_0, -1, Nei, Index_Nei);
  printf(" (1) There are %d neighboring configurations around ", a);
  Writing_General_Configuration_Table(D, 1, Configuration);  
  printf(" transitioning UP into it:\n");
  Writing_General_Configuration_Table(D, a, Nei);

  printf(" (2) There are %d neighboring configurations around ", a);
  Writing_General_Configuration_Table(D, 1, Configuration);  
  printf(" transitioning UP into it:\n\n");
  for(i=0; i<a; i++){
    Index = Index_Nei[i];
    Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
  }

  printf(" Defining the network of all configuration states... ... \n");
  getchar();
  /* Defining the network of all configuration states */
  configuration ** Co = (configuration **)calloc(S, sizeof(configuration *));
  for (i=0; i<S; i++) {
    Co[i] = (configuration *)calloc(1, sizeof(configuration));
    Co[i]->n = Configuration_Table[i]; 
    Co[i]->anDw = (int *)calloc(D + 1, sizeof(int));
    Co[i]->anUp = (int *)calloc(D + 1, sizeof(int));   
  }
  Generating_Network_of_Configurations(Co, Configuration_Table, D, S, n_0);
  // Printing_Network_of_Configurations(Co, D, S);

  for (i=0; i<S; i++) {
    free(Co[i]->anDw);
    free(Co[i]->anUp); 
    free(Co[i]);  
  }
  free (Co);
 
  free(Index_Nei); 
  for (i = 0; i < D; i++) 
    free(Nei[i]);
  free(Nei);
 
  free(Configuration[0]);
  free(Configuration);

  for (i = 0; i < S; i++)
    free(Configuration_Table[i]);
  free(Configuration_Table);

  for (i = 0; i < S; i++)
    free(Configuration_Table_D5[i]);
  free(Configuration_Table_D5);

  return 0;
}

void Writing_Configuration_Table(unsigned int D, unsigned int n_0, int ** Configuration_Table)
{
  /* Counting (brut force) number of configuations for D = 5, and writing them out on std output 
  */
  unsigned long long int S;
  int i, j, k, l, m;
  
  assert(D == 5);
 
  S=0; 
  for(i=0; i<=n_0; i++)
    for(j=0; j<=n_0-i; j++)
      for(k=0; k<=n_0-i-j; k++)
	      for(l=0; l<=n_0-i-j-k; l++)
	        for(m=0; m<=n_0-i-j-k-l; m++){
            Configuration_Table[S][0] = m;
            Configuration_Table[S][1] = l;
            Configuration_Table[S][2] = k;
            Configuration_Table[S][3] = j;
            Configuration_Table[S][4] = i;

	          printf("Configuration[%lld] = [%d, %d, %d, %d, %d]\n",
		        S, m, l, k, j, i);
	      
	          S++;  
          }     
} 

unsigned long long int Func_Hyper_Brut_Force_D2(unsigned int D,
						                                    unsigned int n_0)
{
  /* Counting bumber of configuations (brute force) for D=2 */
  unsigned long long int S;
  int i, j;
  
  assert(D == 2);
  
  S=0; 
  for(i=0; i<=n_0; i++)
    for(j=0; j<=n_0-i; j++)
      S++;

  return(S);
}

unsigned long long int Func_Hyper_Brut_Force_D5(unsigned int D,
						                                    unsigned int n_0)
{
  /* Counting bumber of configuations (brute force) for D=5 */
  unsigned long long int S;
  int i, j, k, l, m;
  
  assert(D == 5);
  
  S=0; 
  for(i=0; i<=n_0; i++)
    for(j=0; j<=n_0-i; j++)
      for(k=0; k<=n_0-i-j; k++)
	      for(l=0; l<=n_0-i-j-k; l++)
	        for(m=0; m<=n_0-i-j-k-l; m++)
	          S++;

  return(S);
}

void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0)
  {
    int Vector_of_Adjacent_Configurations(int * , int, int , int , int ** , int *);
    int i, j, a, Index; 

    int ** Nei = (int **)calloc(D, sizeof(int *));
    for(i=0; i<D; i++)
      Nei[i] = (int *)calloc(D, sizeof(int));

    for (i=0; i<S; i++) {
      Co[i]->n = Configuration_Table[i];

      a = Vector_of_Adjacent_Configurations(Co[i]->n, D, n_0, +1, Nei, Co[i]->anUp);
      Co[i]->anUp[D] = a; 

      printf(" There are %d neighboring configurations around c[%d]=[ ", a, i);
      for(j=0; j<D; j++) printf("%d ", Co[i]->n[j]);
      printf("] transitioning down into it:\n");
      for(j=0; j<a; j++){
        Index = Co[i]->anUp[j];
        Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
      }

      a = Vector_of_Adjacent_Configurations(Co[i]->n, D, n_0, -1, Nei, Co[i]->anDw);
      Co[i]->anDw[D] = a; 

      printf(" There are %d neighboring configurations around c[%d]=[ ", a, i);
      for(j=0; j<D; j++) printf("%d ", Co[i]->n[j]);
      printf("] transitioning up into it:\n");
      for(j=0; j<a; j++){
        Index = Co[i]->anDw[j];
        Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
      } 

      printf("...\n\n");
    }

    for(i=0; i<D; i++)
      free(Nei[i]);
    free(Nei);
  }