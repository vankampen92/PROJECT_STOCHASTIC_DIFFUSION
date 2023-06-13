#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct treenode
{
    /* In this implementation, only true tree leaves are ordered      */
    /* The order of internal nodes is not initialized                 */
    int    level; /* Tree Level (0: root; n: leaves)                  */
    int    order; /* Order within the leaf level, from 0 to 2^n - 1   */
    double value;
     
    struct treenode * left;
    struct treenode * right;
    struct treenode * parent; 
}treenode;

treenode * createtreenode (double value, treenode * parent, int level)
{
    treenode * result = (treenode *)malloc(sizeof(treenode));

    if( result != NULL ) {

        result->value  = value; 
        result->level  = level;     /* Result is created at a level */

        result->left   = NULL;
        result->right  = NULL;
        
        result->parent = parent;
    }

    return result;
}

void printtreenode(treenode * Node)
{
    printf("Level: %d\t", Node->level);
    printf("Value: %g\n", Node->value);
}

void printtabs (int No_of_TABS)
{
    int i; 

    for(i=0; i<No_of_TABS; i++) 
        printf("\t");
}

void printtree_rec (treenode * root, int level)
{
    if (root == NULL) {
        printtabs (level);
        printf("--<empty>--\n");
        return;
    }
    else {
        if( root->parent != NULL) {
            printtabs(level);
            printf("Level (parent)= %d\n", root->parent->level);
            printtabs(level);
            printf("Value (parent)= %g\n", root->parent->value);
        }
    }
    
    printtabs(level);
    printf("Level (node)= %d\n", root->level);
    printtabs(level);
    printf("Value (node)= %g\n", root->value);
    
    printtabs(level+1);
    printf("left\n");
    printtree_rec(root->left, root->level+1);
    
    printtabs(level+1); 
    printf("right\n");
    printtree_rec(root->right, root->level+1);
    printtabs(level);

    printf("done\n");
    return;  
}

void printtree(treenode * root)
{
    printf("The whole tree from root node:\n");

    printtree_rec(root, 0);
    return;
}

void deleteTree(treenode *root)
{
    if(root == NULL) return;

    /* Delete Left sub-tree */
    deleteTree(root->left);
    /* Delete right sub-tree */
    deleteTree(root->right);
     
    /* At last, delete root node */
    printf("Deleting Node (level = %d): %g\n", root->level, root->value);
    free(root);
     
    return;
}

// Print only leaf nodes from binary tree.
treenode * leafPrint(treenode * root)
{
    if (root == NULL) return NULL;

    if (root->left == NULL && root->right == NULL) {
        printtreenode(root);
        return root;
    }
 
    // Else recursively print left and right
    // subtrees.
    root->left = leafPrint(root->left);
    root->right = leafPrint(root->right);
 
    return root;
}

// Delete only leaf nodes from binary tree.
treenode * leafDelete(treenode * root)
{
    if (root == NULL) return NULL;

    if (root->left == NULL && root->right == NULL) {
        free(root);
        return NULL;
    }
 
    // Else recursively delete in left and right
    // subtrees.
    root->left = leafDelete(root->left);
    root->right = leafDelete(root->right);
 
    return root;
}

// Function for inorder traversal in a BST.
void printinorder(treenode * root)
{
    if (root != NULL) {
        printinorder(root->left);
        printf("%g ", root->value);
        printinorder(root->right);
    }
}

int power_int(int a, int n)
{
  /* An integer 'a' raised to the power of an integer 'n' */
  int i;
  int R;

  R = 1;
  for(i=0; i<n; i++)
    R *= a;

  return R;
}

treenode * createBinaryTree_DiscreteDistribution(treenode ** Leaves, int n)
{
    /* Create 2^{n-1} parents, recursively */
    double S; 
    int i, M; 

    M = power_int(2, n-1);
    treenode ** Parents = (treenode **)malloc(M * sizeof(treenode *));
    for(i=0; i<M; i++) {
        printf(" Parent: %d\n", i);
        S = Leaves[2*i]->value + Leaves[2*i+1]->value;         
        Parents[i] = createtreenode(S, NULL, n-1);
        Parents[i]->left   = Leaves[2*i];
        Parents[i]->right  = Leaves[2*i+1];
        Leaves[2*i]->parent   = Parents[i]; 
        Leaves[2*i+1]->parent = Parents[i];

        printtreenode(Parents[i]);
    } 
    treenode * root = Parents[0]; 

    if( n > 1 ) 
        root = createBinaryTree_DiscreteDistribution(Parents, n-1);
    else  
        return root;  
}

void sum_Delta_upto_Root(treenode * root, treenode * Leaf, double Delta)
{
    Leaf->value += Delta; 

    if(Leaf == root) 
        return;

    Leaf = Leaf->parent; 
    
    sum_Delta_upto_Root(root, Leaf, Delta);
}

int choose_Individual_Event(treenode * root, double x)
{
    treenode * node = root; 

    if(node->left == NULL && node->right == NULL)  
        return node->order;  

    if( x < node->left->value ) {
        choose_Individual_Event(node->left, x);
    }
    else {
        x -= node->left->value;
        choose_Individual_Event(node->right, x);
    }   
}

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

    deleteTree(n1); // n1 is the root of all evil. 

    /* Third: Creating a binary tree to sample and
              update a discrete dirtribution:
              
              f = {f_1, ..., f_M}
    */
    printf("\n");
    printf("Creating a binary tree to sample a discrete probability distribution...\n");
    getchar();

    printf("Creating, first, a discrete distribution of 2^n values:\n");
    printf("Enter n, No of TREE LEVELS, from level 0, the root, to level n, the leaves... ");
    scanf("%d", &n); 

    /* The  drand48() returns nonnegative double-precision floating-point values 
       uniformly distributed over the interval [0.0, 1.0)
    */
    M = power_int(2, n); 
    double * Rates     = (double *)calloc(M, sizeof(double));
    treenode ** Leaves = (treenode **)malloc(M * sizeof(treenode *));
    for (i=0; i<M; i++) {
        Rates[i]  = drand48();
        Leaves[i] = createtreenode(Rates[i], NULL, n);
        Leaves[i]->order = i;      
    }

    for (i=0; i<M; i++) 
        printf( "Rate[%d] = %g\n", i, Leaves[i]->value );

    treenode * root = createBinaryTree_DiscreteDistribution(Leaves, n);       

    printtree(root);

    for (i=0; i<M; i++) 
        printf( "Rate[%d] = %g\n", i, Leaves[i]->value );

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

    printf(" Enter a number to exit (the program will exit)... ... ...");
    scanf("%d", &n); 

    deleteTree(root);

    free(Rates);
    return(0);
}