/* David Alonso's implementation of the spanning-cluster algorithm for cluster 
   labeling and calculation of the distribution of cluster sizes in general networks.

   This code is inspired by Tobin Fricke's code:
   http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html 

   Copyright (c) October 10, 2023, by David Alonso <dalonso@ceab.csicx.es>
   Distributed under the terms of the GNU Public License.

   Modified ... 

   This program is written in the 1999 standard of the C language (C99).  Older C
   compilers will refuse to compile it. You can use a C++ compiler, a C99 compiler,
   or you can modify this code to comply with a previous version of the C standard.
   The GCC compiler supports C99 as of gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
   
   Compile this program with (-g for debugging!!!)

   gcc -Wall -g main.c hk.c node.c -o cluster

*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "hk.h"   /* where the Hoshen-Kopelman algorithm is implemented (hk.c)        */
#include "node.h" /* where the "class" node for generic networks is definded (node.c) */

/* The sample program reads in a matrix from file or standard input (or creates it
 * with a particular occupancy probability), and runs two clustering labelling 
 * algorithms (Hoshen-Kopelman and the spanning-cluster) for comparison.
 *  
 * In addition, the distribution of cluster sizes is also calculated.
 *  
 * The form of the input from standard input is two integers giving the
 * dimensions of the matrix, followed by the matrix elements (with data separated by
 * whitespace). For instance, a sample input file is the following:

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
        clusters (connected components) of occupied/activated nodes. 
     3. Apply Tobin's function: 
                 int clusters = hoshen_kopelman(matrix,m,n);
        to compare results 
     4. Transform the collection of clusters or sub-network components of 
        different sizes back into a matrix representation where clusters are 
        left labeled in increasing order from 1 to the total number clusters
        (cluser_id).
     5. Give as an output the size of every cluster in the population 
        of connected subnetworks of occupied/activated nodes and the 
        corresponding distribution of cluster sizes.
*/
/* Further experiment: the relationship between site occupation probability and 
   the resulting number of clusters (if number of rows and columns are both larger 
   than 10) can be explored. 
*/

int main(int argc, char **argv) 
{
  int i, j, k, n, m, i_x, j_y;
  int ** matrix;
  double p; 

  /* Read in the matrix from standard input or let the program generate a random matrix.

     The whitespace-deliminated matrix input is preceeded by the number of rows and 
     number of columns 
  */
  printf("\n");
  printf(" Cluster labelling in generic networks inspired in Tobin Fricke's code for\n"); 
  printf(" the Hoshen-Kopelman algorithm for cluster labeling\n");
  printf(" Here the spanning-cluster algorithm is tested and compared to Tobin's output\n");
  printf(" on regular squared networks. Notice that Tobin's code does not consider\n");
  printf(" periodic boundary conditions. Here, instead, we prescribe these conditions.\n");
  printf(" Results from the two algorithms will only exactly match when there are no clusters\n");
  printf(" wrapping around the boundaries of the squared grid. For the two algorithms to yield\n");
  printf(" the same result, the boundaries of the square grid should be zeroed.\n");
  printf(" If not, there may be clusters wrapping around the borders, which\n");
  printf(" Hoshen-Kopelman will break them into two different clusters,  \n");
  printf(" while the spanning-cluster algorithm will consider them the same cluster.\n\n");
 
  printf("Enter n (No of Rows)... ");
  scanf("%d", &n); 
  printf("Enter n (No of Columns)... ");
  scanf("%d", &m); 
  // n = rows, m = columns

  printf(" If number of rows and columns are both larger than 10,\n");
  printf(" a random matrix will be generated with an occupancy probability, p\n");
  if( n > 10 && m > 10) {
    printf("Enter p (occupancy probability)... ");
    scanf("%lf", &p);

    matrix = (int **)calloc(n, sizeof(int*));
    for (i=0; i<n; i++) 
      matrix[i] = (int *)calloc(m, sizeof(int));

    for (i=0; i<n; i++) 
      for (j=0; j<m; j++)
	      matrix[i][j] = (drand48() < p);
  }
  else {
    printf("Otherwise, enter the matrix row per row...\n");
    matrix = (int **)calloc(n, sizeof(int*));
    for (i=0; i<n; i++) {
      matrix[i] = (int *)calloc(m, sizeof(int));
      for (j=0; j<m; j++)
	      scanf("%d",&(matrix[i][j]));
    }
  }  

  printf(" --input matrix-- \n");
  print_matrix(matrix, n,m);

  printf(" Press any key to start.\n");
  getchar();
    
  int No_of_NODES = m * n; 
  int No_of_NEIGHBORS = 4;    /* Von Newman Squared Network */
  double X_DIMENSION  = 1.0; 
  double Y_DIMENSION  = 1.0;

  node ** nw = (node **)malloc( No_of_NODES * sizeof(node *) );
  Network_Allocation(nw, No_of_NODES, No_of_NEIGHBORS); 

  /* Process the matrix into a proper network of connected nodes */
  initiating_squared_grid_network_from_matrix(matrix, n, m, X_DIMENSION, Y_DIMENSION, 
                                              No_of_NEIGHBORS, nw);
  int max_cluster_size;                                               
  int cluster_size; 
  int * cluster_size_distribution = (int *)calloc(No_of_NODES, sizeof(int));
  /* If the whole network is activated (nw[i]->val = 1 for all i), then there
     will a single cluster of size No_of_NODES: 
                     cluster_size_distribution[No_of_NODES-1] = 1
                     cluster_size_distribution[i] = 0 for all i<(No_of_NODES-1)
                     represening all clusters of sizes smaller than No_of_NODES;
     Here we prescribe: 
          cluster_size_distribution[0]: No of clusters of size 1
          ...
          cluster_size_distribution[k-1]: No of clusters of size k
          ...
          cluster_size_distribution[No_of_NODES-1]: No of clusters of size No_of_NODES            
  */
  /* D will store pointers to the whole population of clusters of different sizes. 
     For instance, { D[1][j][0], D[1][j][1] } would represent a set of two pointers  
     to the 2 nodes corresponding to the j-th cluster of size 2. These pointers point
     to two connected and activated nodes, say, { nw[3], nw[4] }, which make up 
     a particular 2-node cluster. 
  */
  node **** D = (node ****)calloc(No_of_NODES, sizeof(node ***));
  for(k=0; k<No_of_NODES; k++){
    D[k] = (node ***)calloc(No_of_NODES/(k+1), sizeof(node **));
    for (j=0; j < No_of_NODES/(k+1); j++)
      D[k][j] = (node **)calloc(k+1, sizeof(node *));
  }

  int No_of_CLUSTERS = 0;
  int acc_csz        = 0;
  max_cluster_size = 0;   
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
    
  printf("Spanning-cluster algorithm reports %d clusters found\n", 
          No_of_CLUSTERS);
  printf("Hoshen-Kopelman algorithm reports %d clusters found\n", clusters);

  /* Checking the structure D pointing to the different clusters of different 
     sizes. Transform the collection of clusters (or sub-network components) of 
     different sizes back again into a matrix representation where clusters are 
     labeled in increasing order from 1 to the total number of clusters.
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

  printf(" The spanning cluster algorithm reports a total of %d clusters across all sizes\n", 
          No_of_CLUSTERS);
  printf(" The maximum cluster spans over %d nodes.\n", max_cluster_size);
  printf(" The distribution of cluster sizes is: [ ");
  for(i = 0; i<max_cluster_size; i++) 
    if(cluster_size_distribution[i] > 0) 
      printf("[n(%d) = %d] ", i+1, cluster_size_distribution[i]);
  printf("]\n");

  free(cluster_size_distribution);

  Network_Free( nw, No_of_NODES ); 

  for(k=0; k<No_of_NODES; k++){
    for (j=0; j < No_of_NODES/(k+1); j++)
      free(D[k][j]);
    free(D[k]); 
  }
  free(D);

  for (i=0; i<n; i++)
    free(matrix[i]);
  free(matrix);

  return (0);
}
