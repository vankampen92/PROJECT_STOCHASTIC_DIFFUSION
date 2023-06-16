#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct treenode
{
    int value;
     
    struct treenode * left;
    struct treenode * right;
}treenode;

treenode * createtreenode (int value)
{
    treenode * result = (treenode *)malloc(sizeof(treenode));

    if( result != NULL ) {
        result->value = value; 
        result->left  = NULL;
        result->right = NULL;
    }

    return result;
}

void printtabs (int No_of_TABS)
{
    int i; 

    for(i=0; i<No_of_TABS; i++) 
        printf("\t");
}

bool insertnumber(treenode ** rootptr, int value)
{
    treenode * root = * rootptr;

    if( root == NULL) {
        /* Tree Empty */
        * rootptr = createtreenode (value);
        return true; 
    }
    if( value == root->value){
        /* Do nothing */
        return false; 
    }
    if( value < root->value){
        return insertnumber( &(root->left), value);
    }
    else {
        return insertnumber( &(root->right), value);
    }
}

bool findnumber (treenode * root, int value)
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

    printtabs(level);
    printf("value = %d\n", root->value);
    printtabs(level);

    printf("left\n");
    printtree_rec(root->left, level+1);
    printtabs(level); 

    printf("right\n");
    printtree_rec(root->right, level+1);
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

void Upper_Closest_Power_of_Two (int No_of_CELLS, int * No_of_LEAVES, int * No_of_TREE_LEVELS) 
{
  int i;

  * No_of_LEAVES = 0; 
  * No_of_TREE_LEVELS = 0; 

  /* Determine the value of No_of_LEAVES and No_of_TREE_LEVELS */
  if (No_of_CELLS > 1) {
    i = 0; 
    while( No_of_CELLS < power_int(2, i) || No_of_CELLS > power_int(2, i+1)) {
      i++;  
    }

    * No_of_LEAVES      = power_int(2, i+1);
    * No_of_TREE_LEVELS = i+1;
  }
}

int main()
{
    int n; 

    /* First: Creating a Binary Tree */
    treenode * n1 = createtreenode(10);
    treenode * n2 = createtreenode(12);  
    treenode * n3 = createtreenode(15);
    treenode * n4 = createtreenode(17);
    treenode * n5 = createtreenode(20);
    treenode * n6 = createtreenode(25);
    treenode * n7 = createtreenode(29);
    
    n1->left = n2; n1->right = n3;
 
    n2->left = n4; n2->right = n5;

    n3->left = n6; n3->right = n7; 

    printtree(n1);

    free(n7);
    free(n6);
    free(n5);
    free(n4);
    free(n3);
    free(n2);
    free(n1);

    printf("Creating an ordered binary tree...\n");
    getchar();
    /* Second: Creating an ordered binary tree */

    treenode * root = NULL;

    insertnumber( &root, 10);
    insertnumber( &root, 17);
    insertnumber( &root, 15);
    insertnumber( &root, 29);
    insertnumber( &root, 20);
    insertnumber( &root, 9);
    insertnumber( &root, 49);

    printtree(root);

    printf("\n%d [%d]\n", 16, findnumber(root, 16));
    printf("\n%d [%d]\n", 29, findnumber(root, 29));


    printf("\n");
    printf("Determining the closest upper power of two of a given integer number: \n");
    getchar();
    printf("Enter n, an integer number, ... ");
    scanf("%d", &n); 

    int No_of_LEAVES = 0; 
    int No_of_TREE_LEVELS = 0; 
    
    Upper_Closest_Power_of_Two (n, &No_of_LEAVES, &No_of_TREE_LEVELS);

    printf(" No of LEAVES = %d\n", No_of_LEAVES);
    printf(" No of TREE LEVELS = %d\n", No_of_TREE_LEVELS);

    return(0);
}