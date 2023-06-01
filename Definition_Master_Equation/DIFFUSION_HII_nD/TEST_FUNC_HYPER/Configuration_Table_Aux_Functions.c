#include <MODEL.h>

void Writing_Configuration_Table(unsigned int D, unsigned int n_0, int ** Configuration_Table);
unsigned long long int Func_Hyper_Brut_Force_D2(unsigned int D, unsigned int n_0);
unsigned long long int Func_Hyper_Brut_Force_D5(unsigned int D, unsigned int n_0);
void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0); 
unsigned long long int Func_Hyper(unsigned int D, unsigned int k);
void Groups(unsigned int D, int D_Max, int * V, int * s);
unsigned long long int General_Configuration_Table(unsigned int D, unsigned int n_0, 
                                                   int ** Configuration_Table);
int Configuration_to_i_Map(int * Configuration, int D, int n_0);
void i_to_Configuration_Map(int ** Configuration_Table, int i, int * Configuration, int D);
void Writing_General_Configuration_Table(int D,
                                         unsigned long long int S, 
                                         int ** Configuration_Table);
void Comparing_Two_Configuration_Tables(int D, unsigned long long int S, 
                                        int ** Configuration_Table_D5, 
                                        int ** Configuration_Table);
void Pickup_a_Configuration_from_Table(int D, int Index, int ** Configuration_Table);
int Vector_of_Adjacent_Configurations(int * Configuration, int D, int n_0, int d,  
                                      int ** Nei, int * Index_Nei);
void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0);

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

/* Mapping of configuration states (n_1,..., n_D) into just one index (i) */ 
int Configuration_to_i_Map(int * Configuration, int D, int n_0)
{
  /* Mapping of configuration states (n_1,..., n_D) into just one index (i) */ 
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
  /* Selecting the i-th configuration state from Configuration Table and storing it in 
     Configuration array (n_1,..., n_D) 
  */ 
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

void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0)
  {
    int i, j, a, Index; 

    int ** Nei = (int **)calloc(D, sizeof(int *));
    for(i=0; i<D; i++)
      Nei[i] = (int *)calloc(D, sizeof(int));

    for (i=0; i<S; i++) {
      Co[i]->n = Configuration_Table[i];

      a = Vector_of_Adjacent_Configurations(Co[i]->n, D, n_0, +1, Nei, Co[i]->anUp);
      Co[i]->anUp[D] = a; 

      printf(" There are %d neighboring configurations around c[%d]=[ ", a, i);
      for(j=0; j<D; j++) printf("%d ", Co[i]->n[j]);   /* Configuration[i] = [n_1, ..., n_S] */
      printf("] transitioning down into it:\n");
      for(j=0; j<a; j++){
        Index = Co[i]->anUp[j];
        Pickup_a_Configuration_from_Table(D, Index, Configuration_Table);
      }

      a = Vector_of_Adjacent_Configurations(Co[i]->n, D, n_0, -1, Nei, Co[i]->anDw);
      Co[i]->anDw[D] = a; 

      printf(" There are %d neighboring configurations around c[%d]=[ ", a, i);
      for(j=0; j<D; j++) printf("%d ", Co[i]->n[j]);   /* Configuration[i] = [n_1, ..., n_S] */
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