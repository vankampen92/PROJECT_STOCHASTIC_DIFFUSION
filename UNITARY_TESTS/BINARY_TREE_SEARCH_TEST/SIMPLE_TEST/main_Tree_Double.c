#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct treenode
{
    int    level; 
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
        result->level  = level;  /* Result is created at a level */ 

        result->left   = NULL;
        result->right  = NULL;

        result->parent = parent;
    }

    return result;
}

void deleteTree(treenode *root)
{
    if(root == NULL)
        return;
    /* Delete Left sub-tree */
    deleteTree(root->left);
    /* Delete right sub-tree */
    deleteTree(root->right);
     
    /* At last, delete root node */
    printf("Deleting Node : %g\n", root->value);
    free(root);
     
    return;
}

// Delete leaf nodes from binary search tree.
treenode * leafDelete(treenode * root)
{
    if (root == NULL) return NULL;

    if (root->left == NULL && root->right == NULL) {

        printf("Deleting Leaf Node (Value = %g)\n", root->value);

        free(root);
        return NULL;
    }
 
    // Else recursively delete in left and right
    // subtrees.
    root->left = leafDelete(root->left);
    root->right = leafDelete(root->right);
 
    return root;
}

void printtabs (int No_of_TABS)
{
    int i; 

    for(i=0; i<No_of_TABS; i++) 
        printf("\t");
}

bool insertnumber(treenode ** rootptr, treenode * parent, double value, int level)
{
    treenode * root = * rootptr;

    if( root == NULL) {
        /* Tree Empty */
        * rootptr = createtreenode (value, parent, level);
        return true; 
    }
    if( value == root->value){
        /* Do nothing */
        return false; 
    }
    if( value < root->value){
        return insertnumber( &(root->left), root, value, root->level+1);
    }
    else {
        return insertnumber( &(root->right), root, value, root->level+1);
    }
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
 
bool findnumber (treenode * root, double value)
{
    if( root == NULL) return false;

    if( root->value == value ) return true;

    if( root->value > value )
        return findnumber(root->left, value);
    else
        return findnumber(root->right, value);
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

int main()
{
    int i; 
    int M; /* Tree Levels (from level 0 or root level to 
                           level M or leaves level)
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
    
    /* Second: Creating an ordered binary tree */
    printf("Creating an ordered binary tree...\n");
    getchar();
    
    treenode * root = NULL;

    insertnumber( &root, NULL, 10.0, 0);
    insertnumber( &root, NULL, 17.0, 0);
    insertnumber( &root, NULL, 15.0, 0);
    insertnumber( &root, NULL, 29.0, 0);
    insertnumber( &root, NULL, 20.0, 0);
    insertnumber( &root, NULL, 35.0, 0);
    insertnumber( &root, NULL, 29.2, 0);

    printtree(root);

    printf("\n%g [%d]\n", 16.0, findnumber(root, 16.0));
    printf("\n%g [%d]\n", 29.0, findnumber(root, 29.0));

    deleteTree(n1); // n1 is the root of all evil. */

    return(0);
}