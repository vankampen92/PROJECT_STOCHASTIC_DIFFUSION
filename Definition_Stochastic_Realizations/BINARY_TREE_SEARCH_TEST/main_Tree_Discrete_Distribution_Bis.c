#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "treenode.h"

/* This function tests first certain basic functions the the 'treenode.c' library and 
   creates a binary tree to sample discrete distributions by independently allocating 
   first all Parent (internal nodes) and Leave (outer) levels and then creating the 
   binary structure between Parents and Leaves at each level 
*/

int main()
{
    int i; 
    int M; /* No of TREE LEAVES */
    int n; /* No of TREE LEVELS 
              (from level 0 or 'root level' to level n or 'leaf level')
           */
    /* First: Creating a Binary Tree */
    treenode * n1 = createtreenode(10.0, NULL, 0);
    treenode * n2 = createtreenode(12.0, n1, 1);  
    treenode * n3 = createtreenode(15.0, n1, 1);
    treenode * n4 = createtreenode(17.0, n2, 2);
    treenode * n5 = createtreenode(20.0, n2, 2);
    treenode * n6 = createtreenode(25.0, n3, 2);
    treenode * n7 = createtreenode(29.0, n3, 2);
    
    n1->left = n2; n1->right = n3;
 
    n2->left = n4; n2->right = n5;

    n3->left = n6; n3->right = n7; 

    printtree(n1);
    
    printf("Print values in order: ");
    printinorder(n1);
    printf("\n");

    deleteTree(n1); // n1 is the tree root (of all evil). 

    /* Second: Creating a binary tree to sample and
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

    printf(" Enter 0 to exit (the program will exit)... ... ...");
    scanf("%d", &n); 

    Binary_Tree_Free ( root, Leaves, Parent, 
                       No_of_CELLS ); 
    
    // deleteTree(root);

    free(Rates);
    return(0);
}