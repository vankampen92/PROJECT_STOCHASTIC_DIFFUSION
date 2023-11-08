/* David Alonso's implementation of the spanning-clustering algorithm for cluster 
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
   
   Compile this program with (-g- for debugging!!!)

   gcc -Wall -g -std=c99 -c node.c

*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "node.h" /* where the "class" node for generic networks is definded */

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
     neighbors that have not been added yet until no neighbor is activated
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
      . csz, Total No of Nodes of the Spanning Cluster that nucleated
        around the n-th input node (cluster size).
      . D[csz-1][j] is the structure pointing to the j-th cluster of 
      size 'csz'. Through the different calls from main() this structure 
      will end up storing the family of all clusters of different sizes. 
  */
  int k, i;  
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
  
  /* Adding the 'k-th' subnetwork component of size 'csz' into the whole 
     population of clusters of different sizes 
  */
  k = csd[csz-1];
  for(i=0; i<csz; i++){
    D[csz-1][k][i]                = sc[i];
    D[csz-1][k][i]->cluster_class = csz; 
    D[csz-1][k][i]->cluster_id    = cluster_id;
  }

  free(sc); 

  return csz;
}


