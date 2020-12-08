#include "GSL_stat.h"

int da_vector_isnan( double *y, int n)
{
  int i, NOT_a_NUMBER;

  NOT_a_NUMBER = 0;
  
  for(i = 0; i < n; i++) {
    if ( gsl_isnan( y[i]) == 1) NOT_a_NUMBER++;
  }

  return(NOT_a_NUMBER);
}

int da_gsl_vector_isnan(const gsl_vector * y)
{
  double x;
  int i, NOT_a_NUMBER;

  NOT_a_NUMBER = 0;
  
  for(i = 0; i < y->size; i++) {
    x = gsl_vector_get(y, i);
    if ( gsl_isnan(x) == 1) NOT_a_NUMBER++;
  }

  return(NOT_a_NUMBER);
}

bool isDisjoint(const int * set1, int n1, const int * set2, int n2)
{
  /* Input:

     Two sets, for instance:

     . set1: (of length n1)
     . set2: (of lenght n2)

     Output:

     . true, if both sets are disjoint, false otherwise
  */
  int i1, i2;

  int * set1_D = (int *)calloc(n1, sizeof(int) );
  int * set2_D = (int *)calloc(n2, sizeof(int) );
  
  for (i1 = 0; i1<n1; i1++) set1_D[i1] = set1[i1];
  for (i2 = 0; i2<n2; i2++) set2_D[i2] = set2[i2];
  
  void * a1 = (void *)set1_D;
  void * a2 = (void *)set2_D;

  /* Ordering the sets in ascending order */
  gsl_heapsort (a1, n1, sizeof(int), compare_integers);
  gsl_heapsort (a2, n2, sizeof(int), compare_integers);

  i1 = 0; i2 = 0;

  while (i1 < n1 && i2 < n2) {
    if      (set1_D[i1] < set2_D[i2])   i1 = i1 + 1;
    else if (set2_D[i2] < set1_D[i1])   i2 = i2 + 1;
    else {
      free(set1_D);
      free(set2_D);
      return (false);
    }
  }

  free(set1_D);
  free(set2_D);
  return (true);
}

int compare_integers (const void * a, const void * b)
{
  int * A = (int *)a;
  int * B = (int *)b;

  if (*A > *B)
    return 1;
  else if (*A < *B)
    return -1;
  else
    return 0;
}
