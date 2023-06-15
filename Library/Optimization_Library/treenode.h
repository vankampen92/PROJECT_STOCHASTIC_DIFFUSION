#ifndef BINARY_TREE_STRUCTURE 
  #define BINARY_TREE_STRUCTURE
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
int choose_Individual_Event(treenode * root, double x);