/* David Alonso's implementation of the spanning clustering algorithm for cluster 
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

   gcc -Wall -g -std=c99 main_SPAGHETTI.c hk.c -o cluster
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
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
  int cluster_id;       /* cluster id label (1, 2, ...)             */
                        /* if cluster_id is 0, the cluster has not 
                           benn labeled yet 
                        */  
  int dg;               /* node degree                              */
  struct node ** AD;
  int * AD_List;        /* List of nodes IDs a node is connected to */
}node;

void  Set_Von_Neumann_1st_Neighbors(node ** PATCH, int N_X, int N_Y, 
                                    int i)
{
  /* This function creates the neighborhood of a node i in an array of nodes 
     and its adjacency list, where nodes are arranged in a N_X times N_Y
     squared grid (with Von Neumann neighborhood and periodoic boundary 
     conditions). 
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
  /* Network allocation: only valid when the number of connections of 
     every node is a constant number, No_of_NEIGHBORS. Alternatively, 
     No_of_NEIGHBORS should be the maximum number of possible neighbors
     any node can have. 
  */
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

void add_node_to_spanning_cluster(node * cn, node ** sc, int * csz, 
                                  int cluster_id)
{
  /* This function looks at the activated nodes in the set of neighbors 
     of a given node 'cn'. If they are activated and not added yet to 
     the spanning cluster, it adds them to the current spanning cluster, 
     'sc' and, recursively, call the same function with each of the activeted 
     neighors that have not been added yet until no neighbor is activated
     (no further adding is possible).
  */
  int i;
  int added_neighbors; 
  
  sc[* csz] = cn;
  cn->cluster_id = cluster_id;
  cn->cluster_class = * csz;  
  (* csz)++;

  added_neighbors = 0; 
  for(i=0; i<cn->dg; i++){
      if (cn->AD[i]->val == 1 && cn->AD[i]->cluster_id == 0) {
        added_neighbors++;
        add_node_to_spanning_cluster(cn->AD[i], sc, csz, cluster_id);
      }
  }

  if(added_neighbors == 0) return;  
}
 
int spanning_cluster( node ** nw, int n, 
                      int No_of_NODES,  
                      int acc_csz, int * csd, 
                      int cluster_id, node **** D ) 
{
  /* Input: 
      . nw, Full Network Structure (set of nodes)
      . n, ID of an activated node (nw[n]->val = 1).
      . No_of_NODES, total No of NODES of the full network
      . acc_csz, accumulated cluster size
      . csd, distribution of cluster sizes: 
             csd[0] number of clusters of size 1
             ...
             csd[k-1] number of clusters of size k
             ...
      . cluster_id, cluster id label in the population of all clusters  
     Ouput:
      . cluster_size, Total No of Nodes of the Spanning Cluster that nucleated
        around the n-th input node.
      . D structure pointing to all the clusters of different sizes. 
  */
  int k, i;  
  int cluster_size;
  int csz;  
  int max_csz;      /* max size of current spanning cluster */  
 
  max_csz = No_of_NODES - acc_csz; 
  
  /* Reserving memmory for the spanning cluster */
  node ** sc = (node **)calloc(max_csz, max_csz); 
  assert(nw[n]->val == 1);
  sc[0] = nw[n];  /* First node in the spanning cluster */
                  /* The spanning cluster grows around this 
                     first node that acts as a 'nucleous' 
                  */
  /* Look at the activated nodes in the set of neighbors of 
     node sc[0]. If they are activated and not added yet to 
     the spanning cluster, add them to the current spanning 
     cluster, 'sc' and, recursively, call the same function with 
     each of the activeted neighors until no neighbor is 
     activated (or, if there still are, they have been already 
     added).
  */
  csz  = 0;  
  add_node_to_spanning_cluster(sc[0], sc, &csz, cluster_id);
  
  /* Storing the 'k-th' subnetwork component of size 'csz' into the whole 
     population of clusters of different sizes 
  */
  k = csd[csz-1];
  for(i=0; i<csz; i++){
    D[csz-1][k][i]  = sc[i];
    D[csz-1][k][i]->cluster_class = csz; 
    D[csz-1][k][i]->cluster_id    = cluster_id;
  }

  free(sc); 

  cluster_size = csz; 
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
     1. Transform the inpupt matrix in a regular squared network or grid
        (Von Neumann neighborhoood), where occupied/activated nodes are 
        labeled with 1, and non-occupied/non-activated nodes with 0.     
     2. Test the spanning cluster algorithm to label and locate all the 
        clusters (connected components) of occupied nodes. 
     3. Transform the collection of clusters or sub-network components of 
        different sizes back into a matrix representation where clusters are 
        left labeled in increasing order from 1 to the total number clusters
        (cluser_id).
     4. Apply Tobin's function: 
                 int clusters = hoshen_kopelman(matrix,m,n);
        to compare results 
     5. Use the output to write the size of every cluster in the population 
        of connected subnetworks of occupied/activaed nodes and the 
        corresponding distribution of cluster sizes.
*/

int main(int argc, char **argv) {

  int i, j, k, n, m, i_x, j_y;
  int ** matrix;

  /* Read in the matrix from standard input

     The whitespace-deliminated matrix input is preceeded
     by the number of rows and number of columns */
  printf("\n");
  printf(" Cluster labelling in generic networks inspired in Tobin Fricke's code for\n"); 
  printf(" the Hoshen-Kopelman algorithm for cluster labeling\n");
  printf(" Here the spanning-cluster algorithm is tested and compared to Tobin's output\n");
  printf(" on regular squared networks. Notice that Tobin's code does not consider\n");
  printf(" periodic boundary conditions. Here, instead, we prescribe these conditions.\n");
  printf(" Results from the two algorithms will only exactly match when there are no clusters\n");
  printf(" wrapping around the boundaries of the squared grid. For the two algorithm to yield\n");
  printf(" the same result, the boundaries of the square grid should be zeroed.\n");
  printf(" If not, there may be clusters wrapping around the borders, which\n");
  printf(" Hoshen-Kopelman will break them into two different clusters,  \n");
  printf(" while the spanning-cluster algorithm will consider them the same cluster.\n\n");
 
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
      for (int j=0; j<m; j++)
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
  int max_cluster_size;                                               
  int cluster_size; 
  int * cluster_size_distribution = (int *)calloc(No_of_NODES, sizeof(int));
  /* If the whole network is activated (nw[i]->val = 1 for all i), then there
     will a single cluster of size No_of_NODES: 
                     cluster_size_distribution[No_of_NODES-1] = 1
                     cluster_size_distribution[i] = 0 for all i<No_of_NODES;
     Here we prescribe: 
          cluster_size_distribution[0]: No of clusters of size 1
          ...
          cluster_size_distribution[k-1]: No of clusters of size k
          ...
          cluster_size_distribution[No_of_NODES-1]: No of clusters of size No_of_NODES            
  */
  /* D will store pointers to the whole population of clusters of 
     different sizes. 
  */
  node **** D = (node ****)calloc(No_of_NODES, sizeof(node ***));
  for(k=0; k<No_of_NODES; k++){
    D[k] = (node ***)calloc(No_of_NODES/(k+1), sizeof(node **));
    for (j=0; j < No_of_NODES/(k+1); j++)
      D[k][j] = (node **)calloc(k+1, sizeof(node *));
  }

  int No_of_CLUSTERS = 0;
  int acc_csz        = 0;  
  for(i=0; i<No_of_NODES; i++) {
    if(nw[i]->val == 1 && nw[i]->cluster_id == 0){
      No_of_CLUSTERS++; /* Increasing labels for cluster_id's 
                           from 1, 2, 3, ... */ 
      cluster_size = spanning_cluster(nw, i, 
                                      No_of_NODES,
                                      acc_csz, cluster_size_distribution, 
                                      No_of_CLUSTERS, D);
      acc_csz     += cluster_size;
      
      max_cluster_size = max(max_cluster_size, cluster_size);

      cluster_size_distribution[cluster_size-1]++; 
    }
  }

  /* Transform the collection of clusters (or sub-network components) of 
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
  printf(" --output (from the Hoshen-Kopelman algorithm)-- \n");    
  print_matrix(matrix,n,m);
    
  printf("Spanning clustering algorithm reports %d spanning clusters found\n", 
          No_of_CLUSTERS);
  printf("Hoshen-Kopelman algorithm reports %d clusters found\n", clusters);

  /* Checking the structure D pointing to the different clusters of different 
     sizes
  */
  for(i = 0; i<max_cluster_size; i++) 
    if(cluster_size_distribution[i] > 0) {
        printf("As many as %d clusters of size %d has been found\n", 
                cluster_size_distribution[i], i+1);
      for(k=0; k<cluster_size_distribution[i]; k++){
        printf("The %d-th cluster of size %d has been processed\n", 
                k+1, i+1);
      for(j=0; j<i+1; j++) {
        i_x = D[i][k][j]->id / n;
	      j_y = D[i][k][j]->id % n;
      
        matrix[i_x][j_y] = D[i][k][j]->cluster_id; 
      }
    }
  }
  /* Output the result */
  printf(" --output (from the spanning cluster algorithm)-- \n");    
  print_matrix(matrix,n,m);

  printf(" The distribution of cluster sizes is: [ ");
  for(i = 0; i<max_cluster_size; i++) 
    if(cluster_size_distribution[i] > 0) 
      printf("[n(%d) = %d] ", i+1, cluster_size_distribution[i]);
  printf("] ");
  
  for(k=0; k<No_of_NODES; k++){
    for (j=0; j < No_of_NODES/(k+1); j++)
      free(D[k][j]);
    free(D[k]); 
  }
  free(D);
  
  for (i=0; i<n; i++)
    free(matrix[i]);
  free(matrix);

  Network_Free( nw, No_of_NODES ); 

  free(cluster_size_distribution);
  
  return 0;
}




