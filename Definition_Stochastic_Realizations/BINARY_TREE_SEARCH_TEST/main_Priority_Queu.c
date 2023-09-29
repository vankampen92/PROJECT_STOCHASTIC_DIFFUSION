#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "treenode.h"

/* This function creates a priority queue using a binariy tree and test the library functions */

int main()
{   
    double l_x; 
    int i, l; 
    int M; /* No of TREE LEAVES */
    int n; /* No of TREE LEVELS 
              (from level 0 or 'root level' to level n or 'leaf level')
           */

    /* First: Creating a priority quee with a binary tree 
    */
    printf("\n");
    printf(" Creating a priority queu, which maintains the minimum at the root level\n");
    printf(" in spite of any change in tree node values. The 'bubbling algorighm' reorders\n");
    printf(" the tree to maintain the correct priority.\n");
    printf(" Press any key to start.\n");
    getchar();

    printf("Enter n (number of events)... ");
    scanf("%d", &n); 

    /* The  drand48() returns nonnegative double-precision floating-point values 
       uniformly distributed over the interval [0.0, 1.0)
    */
    M = n; 
    double * Rates    = (double *)calloc(M, sizeof(double));
    treenode ** Priority = (treenode **)calloc(M, sizeof(treenode *));
    
    int No_of_CELLS = M;
    int No_of_LEAVES;
    int No_of_TREE_LEVELS; 

    No_of_TREE_LEVELS = Calculating_No_of_TREE_LEVELS(No_of_CELLS);
    No_of_LEAVES      = power_int(2, No_of_TREE_LEVELS);
    /* No_of_TREE_LEVELS corresponds to the No of (internal) TREE LEVELS */
    
    treenode ** Leaves; 
    treenode *** Parent; 
    treenode * root = Binary_Tree_Allocation( No_of_LEAVES, &Leaves, &Parent );

    assert(root == Parent[0][0]);

    for (i=0; i<M; i++) { 
        Rates[i]  = -1.0/2.0 * log(drand48());
        printf("Rate[%d] = %g\n", i, Rates[i]);
    }

    Rates[5] = +INFINITY; 

    /* Inserting the Rates in the binary tree and creating the associated 
       priority vector of pointers to tree nodes 
    */
    for (i=0; i<M; i++) 
        Priority_Queu_Insert_Value(i, Rates[i], No_of_TREE_LEVELS, 
                                   Priority, Parent, Leaves);    
    printtree(root);
    leafPrint(root); 
    printf("This tree has been ordered according to a priority queu...\n"); 
    getchar();
    
    double Delta = 14.9;
    treenode * Leaf = Leaves[0]; 
    Leaf->value -= Delta; 

    printtree(root);
    leafPrint(root); 
    printf("The tree has changed, and it is not ordered yet\n"); 
    getchar();
   
    bubbling(Leaf, Priority);

    printtree(root);
    leafPrint(root); 
    printf("Tree has been reorded through the recursive bubbling algorithm...\n"); 
    getchar();
    
    assert(2<M);
    treenode * Node = Priority[2];  
    Node->value -= Delta; 

    printtree(root);
    leafPrint(root); 
    printf("The tree has changed again, and it does not look ordered yet\n"); 
    getchar();
   
    bubbling(Node, Priority);

    printtree(root);
    leafPrint(root); 
    printf("Tree has been reorded again through the recursive bubbling algorithm...\n"); 
    getchar();

    for (i=0; i<M; i++) { 
        printf("Rate[%d] = %g\n", i, Rates[i]);
        printf("Priority[%d] = %g\n", i, Priority[i]->value);
    }

    printf(" Enter 0 to exit (the program will exit)... ... ...");
    scanf("%d", &n); 

    Binary_Tree_Free ( root, Leaves, Parent, No_of_LEAVES ); 
    
    free(Rates);
    free(Priority);

    return(0);
}