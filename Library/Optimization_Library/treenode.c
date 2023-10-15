#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "treenode.h"

treenode * createtreenode (double value, treenode * parent, int level)
{
    treenode * result = (treenode *)malloc(sizeof(treenode));

    if( result != NULL ) {

        result->value  = value; 
        result->level  = level; /* Result is created at a level                */
        result->index  = 0;     /* Index, only when functioning as a Priority Queu 
                                   (coupled to an array of pointers) 
                                */ 
        result->left   = NULL;
        result->right  = NULL;
        
        result->parent = parent;
    }

    return result;
}

void printtreenode(treenode * Node)
{
    printf("Level: %d\t", Node->level);
    printf("Value: %g\t", Node->value);
    printf("Index: %d\n", Node->index);
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
            printtabs(level);
            printf("Index (parent)= %d\n", root->parent->index);
        }
    }
    
    printtabs(level);
    printf("Level (node)= %d\n", root->level);
    printtabs(level);
    printf("Value (node)= %g\n", root->value);
    printtabs(level);
    printf("Index (node)= %d\n", root->index);
    
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
       maintain the discrete distribution ready to be sampled. 
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

void partial_sums_upto_p(treenode * root, int p, int  n, double * S)
{
    /* If the binary tree maintains partial sums to sample 
       a discrete probability distribution, P[0], ..., P[2^N - 1],
       where N = n-1, then this recursive function can be used to 
       sum up to a certain p:

                            S = P[0] + ... + P[p-1]  
       
       Notice that the discrete distribution takes values from 
       p[0] to p[2^N - 1] while the input argument is n = N-1!!! 
    */   
    int a;
    treenode * node = root; 

    if(node->left == NULL && node->right == NULL)  
        return;  

    a = p/power_int(2, n);
    p = p%power_int(2 , n);

    if( a  == 1 ) {
        * S += node->left->value;
        n--; 
        partial_sums_upto_p(node->right, p, n, S);
    }
    else {
        n--;
        partial_sums_upto_p(node->left, p, n, S);
    }   
}

treenode * sumBinaryTree_DiscreteDistribution(treenode *** Parent, 
                                              treenode ** Leaves, int n)
{
    /*  Set up the next parent level (n-1) from leaves (at level n)
        The function recursively sets up the partial sums that will 
        maintain the discrete distribution ready to be sampled. 
    */
    double S; 
    int i, M; 
    treenode * root; 
    treenode ** Parents;

    if( n >= 1) {
        M = power_int(2, n-1);

        Parents = Parent[n-1];
        
        for(i=0; i<M; i++) {
            // printf(" Parent: %d\n", i);
        
            S = Leaves[2*i]->value + Leaves[2*i+1]->value;  
        
            Parents[i]->value  = S; 
            Parents[i]->level  = n-1;     /* Result is created at a level */

            // printtreenode(Parents[i]);
        } 

        root = Parents[0]; 

        if( n > 1 ) {
            root = sumBinaryTree_DiscreteDistribution(Parent, Parents, n-1);
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

treenode * Binary_Tree_Setting_Structure(treenode **** Parent, 
                                         treenode *** Leaves, int n)
{
    /*  Set up the next parent level (n-1) from leaves (at level n) 
        without any further settings neither:
        1. Partial sums that will maintain the discrete distribution ready 
           to be sampled (Gillespie Direct Method).
        nor: 
        2. Priority Queu settings that will maitain a minimum value at the 
           root level (Gillespie Next Reaction Method).

        The process is performed recursively unitl root is reached.    
    */ 
    int i, M; 
    treenode * root; 
    treenode ** Parents;

    if( n >= 1) {
        M = power_int(2, n-1);

        Parents = (*Parent)[n-1];
        
        for(i=0; i<M; i++) {
            // printf(" Parent: %d\n", i);

            Parents[i]->left   = (*Leaves)[2*i];
            Parents[i]->right  = (*Leaves)[2*i+1];
            (*Leaves)[2*i]->parent   = Parents[i]; 
            (*Leaves)[2*i+1]->parent = Parents[i];
        
            // printtreenode(Parents[i]);
        } 

        root = Parents[0]; 

        if( n > 1 ) {
            root = Binary_Tree_Setting_Structure(Parent, &Parents, n-1);
        }
        else {
            return root;
        }
    }
    else {
        root = (*Leaves)[0];
        return root; 
    } 
}

treenode *  Binary_Tree_Allocation (int No_of_CELLS, 
                                    treenode *** Leaves, treenode **** Parent)
{
  /* This function allocates the tree (from root). It creates space for the leaves 
     (the outer and final tree level) and all internal nodes for each tree level 
     (the parents). 

     Then, it creates the binary structure of the tree by calling the function: 
    
     Binary_Tree_Setting_Structure();

     which connects parents and leaves, as required, recursively. 

     Output: 

     . root, a pointer to the tree root.  
  */
  int i, k, No_of_LEAVES, No_of_TREE_LEVELS, No;

  /* Determine the value of No_of_LEAVES and No_of_TREE_LEVELS */
  No_of_TREE_LEVELS = 0;   /* Only the root!!! */
  No_of_LEAVES      = 1;   /* The root!!!      */
  
  i = 0; 
  if (No_of_CELLS > 1) {
    while( No_of_CELLS < power_int(2, i) || No_of_CELLS > power_int(2, i+1)) {
      i++;  
    }
    No_of_LEAVES      = power_int(2, i+1);
    No_of_TREE_LEVELS = i+1;
  }

  /* No_of_TREE_LEVELS is the No of (internal) TREE LEVELS 
     (without counting the final outer LEAVE level)
     This number will be used to create all Parent tree levels 
  */

  (*Leaves) = (treenode **)malloc(No_of_LEAVES * sizeof(treenode *));
  for(i=0; i<No_of_LEAVES; i++){ 
      (*Leaves)[i] = createtreenode(0.0, NULL, No_of_TREE_LEVELS);
      (*Leaves)[i]->order = i;
  } 

  (*Parent) = (treenode ***)malloc(No_of_TREE_LEVELS * sizeof(treenode **));
  for(k=0; k<No_of_TREE_LEVELS; k++){ 
    No      = power_int(2, k);  /* Number of Leaves at level k */
    (*Parent)[k] = (treenode **)malloc(No * sizeof(treenode *));
    for(i=0; i<No; i++){ 
      (*Parent)[k][i] = createtreenode(0.0, NULL, k);
      (*Parent)[k][i]->order = i;
    }
  }

  treenode * root = Binary_Tree_Setting_Structure(Parent, Leaves, No_of_TREE_LEVELS);

  return(root);
}

void Binary_Tree_Free ( treenode * root, treenode ** Leaves, treenode *** Parent, 
                        int No_of_CELLS ) 
{
  int i, k, No_of_LEAVES, No_of_TREE_LEVELS; 

  /* Determine the value of No_of_LEAVES and No_of_TREE_LEVELS */
  No_of_TREE_LEVELS = 0;   /* Only the root!!! */
  No_of_LEAVES      = 1;   /* The root!!!      */
  
  i = 0; 
  if (No_of_CELLS > 1) {
    while( No_of_CELLS < power_int(2, i) || No_of_CELLS > power_int(2, i+1)) {
      i++;  
    }
    No_of_LEAVES      = power_int(2, i+1);
    No_of_TREE_LEVELS = i+1;
  }

  deleteTree(root);
  free(Leaves);  

  for(k=0; k<No_of_TREE_LEVELS; k++){ 
    free(Parent[k]);
  }
  free(Parent); 
}

int Calculating_No_of_TREE_LEVELS(int M)
{   
    /* This function returns the label, 'n", of the tree level that correspond to 'M', 
       where 'M' is a node label (0, 1, ...). The output for this function should be
       an integer between 0 and the labe of tree leaves                               */
    /* In fact, the total number of levels in the the tree is always the label of tree
       leaves plus 1, because we have labeled the root as level 0 !!!

        Input: 
        . M, is an input index from equal to 0 to equal to 2*(2^n -1). 
        
       The algorithm determins 'n' from the input value 'M'. 
    */
    int No_of_TREE_LEVELS; 
    int i; 

    if (M == 0) 
        No_of_TREE_LEVELS = 0; /* Only the root level */
    else {
        i = 0; 
        while( M < ( power_int(2, i)-1 ) || M > ( 2*(power_int(2, i)-1) ) ) {
                i++;  
        }
        No_of_TREE_LEVELS = i;
    }
    return(No_of_TREE_LEVELS);
}

void Priority_Queu_Insert_Value( int i, double Rate, int LEAVES_LEVEL, 
                                 treenode ** Priority, 
                                 treenode *** Parent, 
                                 treenode ** Leaves)
{    
    /* This function inserts value 'Rate' corresponding to an index 'i' from 
       i=0 to i=(2^{LEAVES_LEVEL} -1) into a priority queu using a binary tree, 
       where Priority is an array of 'node' pointers, where, specifically, 
       Priority[i] points to the tree node that contains index 'i'. 
       
       The re-use of this Priority Queu without re-allocating the binary tree 
       requires that consecutive nodes should be inserted at increasing values 
       of the input argument, the index 'i'. In this way, the whole tree is     
       correctly re-initialized from the root node up to the leaves with brand 
       new values. 
    */
    int k, n, Sn;     
    /* n, tree level from 0 to LEAVES_LEVEL */ 

    n = Calculating_No_of_TREE_LEVELS(i); /* Tree level from 0 to l-1,  where 
                                             l correspond to the total number
                                             of tree levels, labeled from 
                                             n = 0 (the root) to n = l-1 (the leaves) 
                                          */
    Sn = power_int(2, n) - 1; 
    k = i - Sn;

    if (n == 0) {
            Parent[0][0]->value = Rate; 
            Parent[0][0]->level  = 0; /* Result is created at a level 0 */
            Parent[0][0]->index  = i;
            Priority[i]          = Parent[0][0]; 
    }
    else if (n > 0 && n < LEAVES_LEVEL ) { 
            Parent[n][k]->value = Rate;
            Parent[n][k]->level = n;
            Parent[n][k]->index = i;
            Priority[i]         = Parent[n][k]; 
            bubbling_up(Parent[n][k], Priority);
    }
    else {
            Leaves[k]->value = Rate;
            Leaves[k]->level = n;
            Leaves[k]->index = i;
            Priority[i]      = Leaves[k];
            bubbling_up(Leaves[k], Priority); 
    }    
}

void bubbling(treenode * Node, treenode ** Priority)
{
    treenode * Child_Min;

    if (Node->parent == NULL){ /* Root level has been reached */
        Child_Min = Determining_Child_Min(Node);
        
        if (Child_Min == NULL) return; /* Leaf level has been reached 
                                          (single node tree !!!)  
                                       */
        if (Node->value > Child_Min->value) {
            swap_Node_values(Node, Child_Min, Priority);
            bubbling(Child_Min, Priority);
        }
    }
    else if(Node->value < Node->parent->value ){
          swap_Node_values(Node, Node->parent, Priority); 
          bubbling(Node->parent, Priority);
    }
    else{
        Child_Min = Determining_Child_Min(Node);
        
        if (Child_Min == NULL) return; /* Leaf level has been reached */ 

        if (Node->value > Child_Min->value) {
            swap_Node_values(Node, Child_Min, Priority);
            bubbling(Child_Min, Priority);
        }
    }
}

treenode * Determining_Child_Min(treenode * Node)
{
    treenode * child; 

    if(Node->left == NULL && Node->right == NULL)
        return(NULL);
        
    if(Node->left->value > Node->right->value)
        child = Node->right; 
    else    
        child = Node->left; 

    return(child); 
}

void bubbling_up (treenode * Node, treenode ** Priority)
{
    if (Node->parent == NULL) return; 
    else {
        if(Node->value < Node->parent->value ){
          swap_Node_values(Node, Node->parent, Priority); 
          bubbling_up(Node->parent, Priority);
        }
    }
}

void swap_Node_values(treenode * Node_0, treenode * Node_1, treenode ** Priority) 
{
    double value; 
    int index; 
    treenode * Node; 

    Node  = Priority[Node_0->index];
    value = Node_0->value;
    index = Node_0->index;
    
    Priority[Node_0->index] = Priority[Node_1->index];
    Node_0->value = Node_1->value;
    Node_0->index = Node_1->index;

    Priority[Node_1->index] = Node;
    Node_1->value = value;
    Node_1->index = index;
}
