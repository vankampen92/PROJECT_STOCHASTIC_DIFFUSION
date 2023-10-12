/* David Alonso's implementation of the Hoshen-Kopelman algorithm for cluster 
   labeling in general networks.

   This code is inspired by Tobin Fricke's code:
   http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html 

   Copyright (c) October 10, 2023, by David Alonso <dalonso@ceab.csicx.es>
   Distributed under the terms of the GNU Public License.

   Modified ... 

   This program is written in the 1999 standard of the C language (C99).  Older C
   compilers will refuse to compile it. You can use a C++ compiler, a C99 compiler,
   or you can modify this code to comply with a previous version of the C standard.
   The GCC compiler supports C99 as of gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
   
   Compile this program with (-g- for debugging!!!)

   gcc -Wall -g -std=c99 main.c hk.c -o cluster

*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
/* #include "node.h" */ /* where the "class" node is definded */
#include "hk.h"

struct point
{
  float x;
  float y;
};

typedef struct node
{
  struct point center;  /* Coordinates of the position of the center of patch */

  int id;               /* node identification                      */
    
  int val;              /* 1 or 0                                   */

  int cluster_class;    /* size of the cluster the node belongs to  */
  int cluster_id;       /* cluster identification                   */  
     
  int dg;               /* node degree                              */
  struct node ** AD;
  int * AD_List;        /* List of nodes IDs a node is connected to */
}node;

void  Set_Von_Neumann_1st_Neighbors(node ** PATCH, int N_X, int N_Y, 
                                    int i)
{
  /* This function creates the neighborhood a node i in a list and its adjacency list, 
     where nodes are initially arranged in a squared matrix N_X times N_Y 
  */
  int i_x, i_y;
  int n_x, n_y;
  int nei;

  i_x = i/N_X;
  i_y = i%N_X;

  /* Upper Neighbor */
  n_x = i_x;
  n_y = (i_y+1)%N_Y;                    /* Periodic Boundary Condition */
  nei = n_x * N_X  + n_y;
  PATCH[i]->AD[0] = PATCH[nei];
  PATCH[i]->AD_List[0] = nei;

  /* Right Neighbor */
  n_x = (i_x+1)%N_X;                    /* Periodic Boundary Condition */
  n_y = i_y;
  nei = n_x * N_X  + n_y;
  PATCH[i]->AD[1] = PATCH[nei];
  PATCH[i]->AD_List[1] = nei;

  /* Lower Neighbor */
  n_x = i_x;
  n_y = (i_y == 0) ? (N_Y-1) : (i_y-1); /* Periodic Boundary Condition */
  nei = n_x * N_X  + n_y;
  PATCH[i]->AD[2] = PATCH[nei];
  PATCH[i]->AD_List[2] = nei;

  /* Left Neighbor */
  n_x = (i_x == 0) ? (N_X-1) : (i_x-1); /* Periodic Boundary Condition */
  n_y = i_y;
  nei = n_x * N_X  + n_y;
  PATCH[i]->AD[3] = PATCH[nei];
  PATCH[i]->AD_List[3] = nei;
}

void Writing_Adjacency_List(node ** PATCH, int No_of_NODES)
{
  int i,j, no;

  no  = No_of_NODES;
  for(i=0; i<no; i++){
      printf("%s %d %s", "Node No", i, "is conntected to [  ");
      for(j=0; j<PATCH[i]->dg; j++) {
        printf("%d ", PATCH[i]->AD_List[j]);
      }
      printf("%s\n", " ].");
    }
    printf("\n");
}

void Writing_Adjacency_List_VonNeumann(node ** PATCH, int No_of_NODES, int N_X, int N_Y)
{
  int i,j, no;
  int i_x, i_y;
  int j_x, j_y;

  no    = No_of_NODES;

  for(i=0; i<no; i++){
    i_x = i/N_X;
    i_y = i%N_X;
    printf("%s %d %s (%d, %d) %s", "Grid node of ID=", i, "located at ", i_x, i_y, "is conntected to [");
    for(j=0; j<PATCH[i]->dg; j++) {
	      j_x = PATCH[i]->AD_List[j]/N_X;
	      j_y = PATCH[i]->AD_List[j]%N_X;
        printf("  [grid node %d located at (%d, %d)]  ",
	              PATCH[i]->AD_List[j], j_x, j_y);
      }
      printf("%s\n", " ].");
    }
    printf("\n");
}

void Network_Allocation(node ** nw, int No_of_NODES, int No_of_NEIGHBORS)
{
  /* Network allocation */
  int i; 

  for(i=0; i<No_of_NODES; i++){
    nw[i] = (node *)calloc(1, sizeof(node));
    nw[i]->AD_List = (int *)calloc(No_of_NEIGHBORS+1, sizeof( int ));
    nw[i]->AD      = (node **)calloc(No_of_NEIGHBORS, sizeof(node *) );
  }
}

void Network_Free( node ** nw, int No_of_NODES)
{
  int i;
      
  for (i=0; i<No_of_NODES; i++) {
    free(nw[i]->AD);
    free(nw[i]->AD_List);
    free(nw[i]);
  }
  free(nw);
}

void initiating_squared_grid_network_from_matrix(int ** matrix, int N_X, int N_Y, 
                                                 double X_DIMENSION, double Y_DIMENSION, 
                                                 int No_of_NEIGHBORS, 
                                                 node ** nw)
{ 
  int i, No_of_NODES; 
  int i_x, j_y; 
  double STEP_X, STEP_Y;

  No_of_NODES = N_X * N_Y;  
  STEP_X      = X_DIMENSION/(double)N_X;
  STEP_Y      = Y_DIMENSION/(double)N_Y;
  
  for(i=0; i<No_of_NODES; i++){

	  i_x = i/N_X;
	  j_y = i%N_X;

    nw[i]->center.x = (double)j_y + 0.5*STEP_X;
	  nw[i]->center.y = (double)i_x + 0.5*STEP_Y;

    if(matrix[i_x][j_y] == 1) 
      nw[i]->val = 1; 
    else
      nw[i]->val = 0; 
	  
    nw[i]->id = i; 
    nw[i]->cluster_class = 0; 
    nw[i]->cluster_id    = 0;
    nw[i]->dg = No_of_NEIGHBORS; 
    nw[i]->AD_List[No_of_NEIGHBORS] = No_of_NEIGHBORS;   

	  Set_Von_Neumann_1st_Neighbors(nw, N_X, N_Y, i);
  }
}

int spanning_cluster(node ** nw, int No_of_NODES, int n) 
{
  /* Input: 
      . nw, Full Network Structure 
      . No_of_Nodes, 
      . n, ID of an activated node (nw[n]->val = 1).
     Ouput:
      . cluster_size, No of Nodes of the Spanning Cluster starting at the n-th node.
  */ 
  int cluster_size; 

  cluster_size = 1; 

  assert(nw[n]->val == 1);

  return cluster_size;
}

/* The sample program reads in a matrix from file or standard input, and runs 
   the clustering labelling algorithm. After that, it gives the distribution of cluster
   sizes 
   
   The form of the input from standard input is two integers giving the
   dimensions of the matrix, followed by the matrix elements (with data separated by
   whitespace). For instance, a sample input file is the following:

   8 8
   1 1 1 1 1 1 1 1
   0 0 0 0 0 0 0 1
   1 0 0 0 0 1 0 1
   1 0 0 1 0 1 0 1
   1 0 0 1 0 1 0 1
   1 0 0 1 1 1 0 1
   1 1 1 1 0 0 0 1
   0 0 0 1 1 1 0 1 

   This sample input gives the following output:

   --input-- 
   1   1   1   1   1   1   1   1 
   0   0   0   0   0   0   0   1 
   1   0   0   0   0   1   0   1 
   1   0   0   1   0   1   0   1 
   1   0   0   1   0   1   0   1 
   1   0   0   1   1   1   0   1 
   1   1   1   1   0   0   0   1 
   0   0   0   1   1   1   0   1 
   --output-- 
   1   1   1   1   1   1   1   1 
   0   0   0   0   0   0   0   1 
   2   0   0   0   0   2   0   1 
   2   0   0   2   0   2   0   1 
   2   0   0   2   0   2   0   1 
   2   0   0   2   2   2   0   1 
   2   2   2   2   0   0   0   1 
   0   0   0   2   2   2   0   1 

   The program reports 2 clusters found of sizes ( 15, 19 )

   Unitary test:
     1. Transform the inpupt matrix in a regular squared network or grid,  
        where occupied nodes are labeled with 1, and non-occupied nodes 
        with 0.     
     2. Test the algorithm to label and locate all the clusters of 
        occupied nodes. 
     3. Transform the collection of clusters or network components of 
        different sizes back into a matrix representation where clusters are 
        labeled in increasing order from 1 to the total number clusters. 
     4. Apply Tobin's function: 
                 int clusters = hoshen_kopelman(matrix,m,n);
        to compare results 
     5. Use the output to write the size of every cluster and the 
          corresponding distribution of cluster sizes.
*/

int main(int argc, char **argv) {

  int i, n, m, i_x, j_y;
  int ** matrix;

  /* Read in the matrix from standard input

     The whitespace-deliminated matrix input is preceeded
     by the number of rows and number of columns */

  printf("\n");
  printf(" Cluster labelling in generic networks based on Tobin Fricke's implementation of\n"); 
  printf(" the Hoshen-Kopelman algorithm for cluster labeling\n");
 
  printf(" Press any key to start.\n");
  getchar();

  printf("Enter n (No of Rows)... ");
  scanf("%d", &n); 
  printf("Enter n (No of Columns)... ");
  scanf("%d", &m); 
  // n = rows, m = columns
  
  printf("Enter the matrix row per row...\n");
  matrix = (int **)calloc(n, sizeof(int*));
  for (int i=0; i<n; i++) {
      matrix[i] = (int *)calloc(m, sizeof(int));
      for (int j=0; j<n; j++)
	      scanf("%d",&(matrix[i][j]));
  }
    
  printf(" --input matrix-- \n");
  print_matrix(matrix, n,m);
    
  int No_of_NODES = m * n; 
  int No_of_NEIGHBORS = 4;    /* Von Newman Squared Network */
  double X_DIMENSION  = 1.0; 
  double Y_DIMENSION  = 1.0;

  node ** nw = (node **)malloc( No_of_NODES * sizeof(node *) );
  Network_Allocation( nw, No_of_NODES, No_of_NEIGHBORS); 

  /* Process the matrix into a proper network of connected nodes */
  initiating_squared_grid_network_from_matrix(matrix, n, m, X_DIMENSION, Y_DIMENSION, 
                                              No_of_NEIGHBORS, nw);
  int cluster_size; 
  int * cluster_size_distribution = (int *)calloc(No_of_NODES+1, sizeof(int));
  /* If the whole network is activated (nw[i]->val = 1 for all i), then there
       will a single cluster of size No_of_NODES: 
                     cluster_size_distribution[No_of_NODES] = 1
                     cluster_size_distribution[i] = 0 for i<No_of_NODES; 
  */
  int No_of_CLUSTERS = 0; 
  for(i=0; i<No_of_NODES; i++) {
    if(nw[i]->val == 1 && nw[i]->cluster_id == 0){
      cluster_size = spanning_cluster(nw, No_of_NODES, i);

      cluster_size_distribution[cluster_size]++; 
      No_of_CLUSTERS++;
    }
  }

  /* Transform the collection of clusters (or network components) of 
     different sizes back into a matrix representation where clusters are 
     labeled in increasing order from 1 to the total number clusters.
  */
  for(i=0; i<No_of_NODES; i++) {
      i_x = i/n;
	    j_y = i%n;

      matrix[i_x][j_y] = nw[i]->val;
    }
  printf(" --check: input matrix (recovered from network)--\n");
  print_matrix(matrix, n,m);

  int clusters = hoshen_kopelman(matrix, n,m);
    /* Output the result */
  printf(" --output-- \n");    
  print_matrix(matrix,n,m);
    
  printf("HK reports %d spanning clusters found\n", No_of_CLUSTERS);
  printf("HK reports %d HK clusters found\n", clusters);

  for (i=0; i<n; i++)
    free(matrix[i]);
  free(matrix);

  Network_Free( nw, No_of_NODES ); 

  free(cluster_size_distribution);
  
  return 0;
}




