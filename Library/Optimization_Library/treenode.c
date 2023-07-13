#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "treenode.h"

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

void deleteTree(treenode * root)
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
    /* Create 2^{n-1} treenode parents from leaves, only if necessary, but 
       always sum leaf level recursively to set up the partial sums that will 
       maintain the discrite distribution ready to be sampled. 
    */
    double S; 
    int i, M; 
    treenode * root; 
    treenode ** Parents;

    if( n >= 1) {
        M = power_int(2, n-1);

        Parents = (treenode **)malloc(M * sizeof(treenode *));
        
        for(i=0; i<M; i++) {
            // printf(" Parent: %d\n", i);
        
            S = Leaves[2*i]->value + Leaves[2*i+1]->value;  
        
            if(Leaves[2*i]->parent == NULL) {       
                Parents[i] = createtreenode(S, NULL, n-1);
                Parents[i]->left   = Leaves[2*i];
                Parents[i]->right  = Leaves[2*i+1];
                Leaves[2*i]->parent   = Parents[i]; 
                Leaves[2*i+1]->parent = Parents[i];
            }
            else {
                Parents[i] = Leaves[2*i]->parent;
                Parents[i]->value = S; 
            } 
            // printtreenode(Parents[i]);
        } 

        root = Parents[0]; 

        if( n > 1 ) {
            root = createBinaryTree_DiscreteDistribution(Parents, n-1);
        }
        else {
            return root;
        }
    }
    else {
        root = Leaves[0];
        return root; 
    } 
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

