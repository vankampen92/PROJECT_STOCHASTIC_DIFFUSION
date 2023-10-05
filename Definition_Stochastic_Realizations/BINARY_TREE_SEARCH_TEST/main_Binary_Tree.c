#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "treenode.h"

/* This function creates a binary tree to sample discrete distributions and 
   test the library functions */

int main()
{
    int i; 
    int M; /* No of TREE LEAVES */
    int n; /* No of TREE LEVELS 
              (from level 0 or 'root level' to level n or 'leaf level')
           */

    /* First: Creating a binary tree to sample and
               update a discrete dirtribution:
              
              f = {f_1, ..., f_M}
    */
    printf("\n");
    printf(" Creating a binary tree to sample a discrete probability distribution...\n");
    printf(" Press any key to start.\n");
    getchar();

    printf("Creating, first, a discrete distribution of 2^n values:\n");
    printf("Enter n, No of TREE LEVELS, from level 0, the root, to level n, the leaves... ");
    scanf("%d", &n); 

    /* The  drand48() returns nonnegative double-precision floating-point values 
       uniformly distributed over the interval [0.0, 1.0)
    */
    M = power_int(2, n); 
    double * Rates     = (double *)calloc(M, sizeof(double));
    
    int No_of_CELLS;
    int No_of_LEAVES;
    int No_of_TREE_LEVELS; 

    No_of_CELLS       = M; 
    No_of_TREE_LEVELS = n;

    treenode ** Leaves; 
    treenode *** Parent; 
    treenode * root = Binary_Tree_Allocation( No_of_CELLS, &Leaves, &Parent );

    for (i=0; i<M; i++) {
        Rates[i]  = drand48();
        Leaves[i]->value = Rates[i];
        Leaves[i]->order = i;      
    }

    root = sumBinaryTree_DiscreteDistribution(Parent, Leaves, No_of_TREE_LEVELS); 

    for (i=0; i<M; i++) 
        printf( "Rate[%d] = %g\n", i, Leaves[i]->value );
      
    printtree(root);
    leafPrint(root); 

    double Delta = 2.0;
    treenode * Leaf = Leaves[2]; 
    sum_Delta_upto_Root(root, Leaf, Delta);
    
    printtree(root);

    for (i=0; i<M; i++) 
        printf( "Rate[%d] = %g\n", i, Leaves[i]->value ); 

    leafPrint(root); 

    double x = root->value * drand48();
    int Event = choose_Individual_Event(root, x);
    printf(" The winner event is... the %d-th!!!\n", Event);
    printf(" (from 0 to %d)\n", M-1);

    double S; 
    int p = (M * drand48());
    printf(" Accummulated distribution up to p = %d\n", p);
    if( p > 0) 
        printf(" S = P[0] + ... + P[%d]\n", p-1);
    else
        printf(" S = 0");

    S = 0.0; 
    partial_sums_upto_p(root, p, n-1, &S);
    printf("S_1 = P[0] + ... + P[%d] = %g\n", p-1, S);

    S = 0.0; 
    for (i=0; i<p; i++) 
        S += Leaves[i]->value; 
    printf("S_2 = P[0] + ... + P[%d] = %g\n", p-1, S);

    printf("S_1 and S_2 should perfectly match!!!\n");

    printf(" If p = 0, then S should be also 0!!!\n");
    S = 0.0;
    partial_sums_upto_p(root, 0, n-1, &S);
    printf("S = %g\n", S);
    if( S == 0.0)
        printf("Your code works perfectly fine!!!\n");

    printf(" Enter 0 to exit (the program will exit)... ... ...");
    scanf("%d", &n); 

    Binary_Tree_Free ( root, Leaves, Parent, 
                       No_of_CELLS ); 
    
    free(Rates);
    return(0);
}