#ifndef BINARY_TREE_STRUCTURE 
  #define BINARY_TREE_STRUCTURE
  typedef struct treenode
  {
    /* In this implementation, only true tree leaves are ordered           */
    /* The order of internal nodes is not initialized                      */
    int    level; /* Tree Level (0: root; n: leaves)                       */
    int    order; /* Order within the leaf level, from 0 to 2^n - 1        */
    int    index; /* Index (when using the binary tree as a priority queu) */
    
    double value;

    struct treenode * left;
    struct treenode * right;
    struct treenode * parent;
  }treenode;
#endif

treenode * createtreenode (double value, treenode * parent, int level);
void printtreenode(treenode * Node);
void printtabs (int No_of_TABS);
void printtree_rec (treenode * root, int level);
void printtree(treenode * root);
void deleteTree(treenode *root);
treenode * leafPrint(treenode * root);
treenode * leafDelete(treenode * root);
void printinorder(treenode * root);
int power_int(int a, int n);
treenode * createBinaryTree_DiscreteDistribution(treenode ** Leaves, int n);
void sum_Delta_upto_Root(treenode * root, treenode * Leaf, double Delta);
void partial_sums_upto_p(treenode * root, int p, int  n, double * S);
int choose_Individual_Event(treenode * root, double x);
treenode * Binary_Tree_Allocation (int No_of_CELLS, 
                                   treenode *** Leaves, treenode **** Parent);
treenode * Binary_Tree_Setting_Structure(treenode **** Parent, 
                                         treenode *** Leaves, int n);
treenode * sumBinaryTree_DiscreteDistribution(treenode *** Parent, 
                                              treenode ** Leaves, int n);
void Binary_Tree_Free ( treenode * root, treenode ** Leaves, treenode *** Parent, 
                        int No_of_CELLS ); 
void Priority_Queu_Insert_Value( int i, double Rate, int LEAVES_LEVEL, 
                                 treenode ** Priority, 
                                 treenode *** Parent, 
                                 treenode ** Leaves);
void bubbling_up (treenode * Node, treenode ** Priority);
void swap_Node_values(treenode * Node_0, treenode * Node_1, treenode ** Priority);
treenode * Determining_Child_Min(treenode * Node);
void bubbling(treenode * Node, treenode ** Priority);
int Calculating_No_of_TREE_LEVELS(int M);


