/* David Alonso's implementation of the spanning-clustering algorithm for cluster 
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
   
   Compile this library with (-g- for debugging!!!)

   gcc -Wall -g -std=c99 -c node.c 
*/
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

void  Set_Von_Neumann_1st_Neighbors(node ** PATCH, int N_X, int N_Y, int i);

void Writing_Adjacency_List(node ** PATCH, int No_of_NODES);

void Writing_Adjacency_List_VonNeumann(node ** PATCH, int No_of_NODES, int N_X, int N_Y);

void Network_Allocation(node ** nw, int No_of_NODES, int No_of_NEIGHBORS);

void Network_Free( node ** nw, int No_of_NODES);

void initiating_squared_grid_network_from_matrix(int ** matrix, int N_X, int N_Y, 
                                                 double X_DIMENSION, double Y_DIMENSION, 
                                                 int No_of_NEIGHBORS, 
                                                 node ** nw);

void add_node_to_spanning_cluster(node * cn, node ** sc, int * csz, int cluster_id);
 
int spanning_cluster( node ** nw, int n, 
                      int No_of_NODES,  
                      int acc_csz, int * csd, 
                      int cluster_id, node **** D ); 

