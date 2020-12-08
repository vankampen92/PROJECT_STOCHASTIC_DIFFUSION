/* Two sets are disjoint set when they have no common elements.
 * In other words, if we get the intersection of two sets, then we will get null set.
 * The method is simple, in this algorithm, two sets are given. The sets should be
 * previously ordered in ascending order. Items are compared between two sets. When
 * there is a match, then it is not a disjoint set, when no items are matched, they
 * are disjoint sets.
 */
#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_heapsort.h>

int compare_integers (const void * a, const void * b);
bool isDisjoint(int * set1, int n1, int * set2, int n2);

/* int main(void) { */
/*   int i; */

/*   int set1[] = {15, 13, 36, 2, 32, 4}; */
/*   int set2[] = {32, 89, 56, 12}; */

/*   printf("Sets to compare before...\n"); */
/*   printf("set1 = { "); */
/*   for(i=0; i<6; i++) printf("%d ", set1[i]); */
/*   printf("}\n"); */

/*   printf("set2= { "); */
/*   for(i=0; i<4; i++) printf("%d ", set2[i]); */
/*   printf("}\n"); */

/*   bool x = true;  /\* equivalent to bool x = 1; *\/ */

/*   x = isDisjoint(set1, 6, set2, 4); */

/*   if (x)  /\* Functionally equivalent to if (x != 0) or if (x != false) *\/ */
/*     puts("set1 and set2 are disjoint!!!\n"); */

/*   if (!x) /\* Functionally equivalent to if (x == 0) or if (x == false) *\/ */
/*     puts("set1 and set2 share an element!!!\n"); */

/*   printf("\nSets to compare after...\n"); */
/*   printf("set1 = { "); */
/*   for(i=0; i<6; i++) printf("%d ", set1[i]); */
/*   printf("}\n"); */

/*   printf("set2 = { "); */
/*   for(i=0; i<4; i++) printf("%d ", set2[i]); */
/*   printf("}\n"); */

/*   return(0); */
/* } */

bool isDisjoint(int * set1, int n1, int * set2, int n2)
{
  /* Input:

     Two sets, for instance:

     . set1: (of length n1)
     . set2: (of lenght n2)

     Output:

     . true, if both sets are disjoint, false otherwise
  */
  int i1, i2;
  void * a1 = (void *)set1;
  void * a2 = (void *)set2;

  /* Ordering the sets in ascending order */
  gsl_heapsort (a1, n1, sizeof(int), compare_integers);
  gsl_heapsort (a2, n2, sizeof(int), compare_integers);

  i1 = 0; i2 = 0;

  while (i1 < n1 && i2 < n2) {
    if (set1[i1] < set2[i2])        i1 = i1 + 1;
    else if (set2[i2] < set1[i1])   i2 = i2 + 1;
    else                         return (false);
  }

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
