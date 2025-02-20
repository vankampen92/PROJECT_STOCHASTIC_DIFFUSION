#include <MODEL.h>

/* This functions allocate, initialize and free a number of local communities,
   (or patches), which make up our total patch system or metacommunity
*/
extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Community_Allocation ( Community ** PATCH, Parameter_Model * P )
{
  int i, j, a;
  int N, no, Sp;

  Sp    = P->LOCAL_STATE_VARIABLES;  /* Total Number of Species 
					                              potentially coexisting locally,
					                              and, therefore, also the
					                              Total Number of State Variables
					                              fully determining the local state
				                             */
  no    = P->No_of_CELLS;

  N     = P->TOTAL_No_of_EVENTS;

  for(i=0; i<no; i++){
    PATCH[i] = (Community *)calloc( 1, sizeof(Community) );

    PATCH[i]->n = (int *)calloc(Sp, sizeof( int ));

    PATCH[i]->rate = (double *)calloc(P->TOTAL_No_of_EVENTS, sizeof( double ));
    PATCH[i]->rToI = (double *)calloc(P->TOTAL_No_of_EVENTS, sizeof( double ));

    PATCH[i]->Patch_Connections = (int *)calloc(P->No_of_NEIGHBORS, sizeof( int ));

    PATCH[i]->NEI  = (Community **)calloc(P->No_of_NEIGHBORS, sizeof(Community *) );

    PATCH[i]->Out_Migration_Vector = (double **)calloc(Sp, sizeof( double *));
    for(a=0; a<Sp; a++)
      PATCH[i]->Out_Migration_Vector[a] = (double *)calloc(P->No_of_NEIGHBORS, sizeof( double ));

    PATCH[i]->In_Migration_Vector = (double **)calloc(Sp, sizeof( double *));
    for(a=0; a<Sp; a++)
      PATCH[i]->In_Migration_Vector[a] = (double *)calloc(P->No_of_NEIGHBORS, sizeof( double ));

    PATCH[i]->Total_Per_Capita_Out_Migration_Rate = (double *)calloc(Sp, sizeof(double) );

    PATCH[i]->Total_Imm_Rate_Preassure = (double *)calloc(Sp, sizeof(double) );

    PATCH[i]->Imm_Rates_Preassure = (double **)calloc(Sp, sizeof( double *));
    for(a=0; a<Sp; a++)
      PATCH[i]->Imm_Rates_Preassure[a] = (double *)calloc(P->No_of_NEIGHBORS, sizeof( double ));

    /* Certain ECOEVO models require the organization of the Event_Delta_Matrix as a tensor,
       where each outter level represents the matrix associated to a species or type 
       of size No_of_EVENTS x No_of_EVENTS, where this is the number of events each species 
       can undergo 
    */
    if(P->TYPE_of_MODEL == 20 || P->TYPE_of_MODEL == 21) {
      PATCH[i]->Event_Delta_Tensor = (double ***)calloc(P->No_of_RESOURCES, sizeof(double **) );
      for(a=0; a<P->No_of_RESOURCES; a++) {
        PATCH[i]->Event_Delta_Tensor[a] = (double **)calloc(P->No_of_EVENTS, sizeof(double *) );
        for(j=0; j<P->No_of_EVENTS; j++)
          PATCH[i]->Event_Delta_Tensor[a][j] = (double *)calloc(P->No_of_EVENTS, sizeof(double) );
      }

      PATCH[i]->Event_Adjacence_List = (int **)calloc(P->No_of_EVENTS, sizeof(int *) );
      for(a=0; a<P->No_of_EVENTS; a++)
        PATCH[i]->Event_Adjacence_List[a] = (int *)calloc(P->No_of_EVENTS+1, sizeof(int) );  
    }
    else {
      PATCH[i]->Event_Delta_Matrix = (double **)calloc(N, sizeof(double *) );
      for(a=0; a<N; a++)
        PATCH[i]->Event_Delta_Matrix[a] = (double *)calloc(N, sizeof(double) );
      
      PATCH[i]->Event_Adjacence_List = (int **)calloc(N, sizeof(int *) );
      for(a=0; a<N; a++)
        PATCH[i]->Event_Adjacence_List[a] = (int *)calloc(N+1, sizeof(int) );
    }

#if defined DIFFUSION_ECO_PLASMIDS
    PATCH[i]->Local_Strain_Population = (Strain **)calloc(P->No_of_RESOURCES, sizeof(Strain *));     
    for(a=0; a<P->No_of_RESOURCES; a++){  
      PATCH[i]->Local_Strain_Population[a] = (Strain *)calloc(1, sizeof(Strain));
    }
    
    PATCH[i]->Local_Plasmid_Population = (Plasmid **)calloc(P->No_of_PLASMIDS, sizeof(Plasmid *));
    for(a=0; a<P->No_of_PLASMIDS; a++)
      PATCH[i]->Local_Plasmid_Population[a] = (Plasmid *)calloc(1, sizeof(Plasmid));

    PATCH[i]->Bacterial_Type_Population = (int *)calloc(P->No_of_STRAINS, sizeof(int)); 
    PATCH[i]->Plasmid_Type_Population   = (int *)calloc(P->No_of_PLASMIDS, sizeof(int));  
#endif
  }
}

void Community_Free (Community ** PATCH, Parameter_Model * P)
{
  int Sp, K, i, j, a;

  Sp  = P->LOCAL_STATE_VARIABLES;     /*  LOCAL_STATE_VARIABLES:
                                          (For instance, 2 times No_of_RESOURCES for
                                          DIFFUSION_ECOEVO_PLANTS model, this is the 
                                          Total Number of species and stages (R or RP, 
                                          for each species) potentially coexisting 
                                          together locally). In general, this is, the
					                                Total Number of State Variables
					                                fully determining the local state
				                             */
  /* BEGIN: Patch Total Destruction */
  for (i=0; i<P->No_of_CELLS; i++){

    free(PATCH[i]->n);

    free(PATCH[i]->rate);
    free(PATCH[i]->rToI);

    free(PATCH[i]->Patch_Connections);

    free(PATCH[i]->NEI);

    for(a=0; a<Sp; a++)
      free(PATCH[i]->Out_Migration_Vector[a]);
    free(PATCH[i]->Out_Migration_Vector);

    for(a=0; a<Sp; a++)
      free(PATCH[i]->In_Migration_Vector[a]);
    free(PATCH[i]->In_Migration_Vector);

    free(PATCH[i]->Total_Per_Capita_Out_Migration_Rate);

    free(PATCH[i]->Total_Imm_Rate_Preassure);

    for(a=0; a<Sp; a++)
      free(PATCH[i]->Imm_Rates_Preassure[a]);
    free(PATCH[i]->Imm_Rates_Preassure);

    if(P->TYPE_of_MODEL == 20 || P->TYPE_of_MODEL == 21) {
      for(a=0; a<P->No_of_RESOURCES; a++) {
        for(j=0; j<P->No_of_EVENTS; j++)
          free(PATCH[i]->Event_Delta_Tensor[a][j]);

        free(PATCH[i]->Event_Delta_Tensor[a]);  
      }
      free(PATCH[i]->Event_Delta_Tensor);

      for(a=0; a<P->No_of_EVENTS; a++)
        free(PATCH[i]->Event_Adjacence_List[a]);

      free(PATCH[i]->Event_Adjacence_List);    
    }
    else {
      for(a=0; a<Sp; a++)                         /*  Sp, LOCAL_STATE_VARIABLES */
        free(PATCH[i]->Event_Delta_Matrix[a]);
      
      free(PATCH[i]->Event_Delta_Matrix);

      for(a=0; a<Sp; a++)                         /*  Sp, LOCAL_STATE_VARIABLES */
        free(PATCH[i]->Event_Adjacence_List[a]);

      free(PATCH[i]->Event_Adjacence_List);
    }  
    
#if defined DIFFUSION_ECO_PLASMIDS
    for(a=0; a<P->No_of_RESOURCES; a++) 
      free(PATCH[i]->Local_Strain_Population[a]);
    free(PATCH[i]->Local_Strain_Population);  

    for(a=0; a<P->No_of_PLASMIDS; a++)
      free(PATCH[i]->Local_Plasmid_Population[a]); 
    free(PATCH[i]->Local_Plasmid_Population);

    free(PATCH[i]->Bacterial_Type_Population);
    free(PATCH[i]->Plasmid_Type_Population); 
#endif

    free(PATCH[i]);
  }

  free( PATCH );

  /*   END: Patch Total Destruction */
}

void Community_Initialization (Community ** PATCH,
			                         Parameter_Model * P )
{
  int i, j, Sp, no;

  Sp  = P->No_of_RESOURCES;
  no  = P->No_of_CELLS;

  for(i=0; i<no; i++){

    PATCH[i]->No_of_RESOURCES = Sp;

    PATCH[i]->No_of_CELLS   = P->No_of_CELLS;
    PATCH[i]->No_of_CELLS_X = P->No_of_CELLS_X;
    PATCH[i]->No_of_CELLS_Y = P->No_of_CELLS_Y;
    PATCH[i]->X_DIMENSION = (double)P->No_of_CELLS_Y;
    PATCH[i]->Y_DIMENSION = (double)P->No_of_CELLS_X;

    PATCH[i]->LOCAL_STATE_VARIABLES  = P->LOCAL_STATE_VARIABLES;
    /* This is the number of dynamic state variables             */
    /* (required to defined the state of the patch)              */
    /* Total Number of Events within a patch: TOTAL_No_of_EVENTS */

    PATCH[i]->ratePatch    = 0.0;      /* Transition probabilities at this patch */
    for( j=0; j<P->TOTAL_No_of_EVENTS; j++) {
      PATCH[i]->rate[j]= 0.0;
      PATCH[i]->rToI[j]= 0.0;
    }

    PATCH[i]->Metapop_Connectivity_Matrix = P->Metapop_Connectivity_Matrix;
  }

  /* When PATCH represents a multi-patch network, patch connections
     and number of patch connections for each patch (PATCH[i]->No_NEI)
     should be initialized                                          */
  Network_Structure_Inititialization (PATCH,
				                              P->No_of_NEIGHBORS,
				                              P->TYPE_of_NETWORK);
#if defined VERBOSE
  Writing_Adjacency_List(PATCH);
  if (P->TYPE_of_NETWORK == 1) Writing_Adjacency_List_VonNeumann(PATCH);
#endif

  Immigration_Preassure_on_Focal_Patch_Initialization( PATCH, P );

#if defined STOCHASTIC_OPTIMIZATION
  /* This optimization can be implemented for any model. So far, 
     it has been implemented for: 
     . MODEL=DIFFUSION_1R1C             TYPE_of_MODEL = 2
     . MODEL=DIFFUSION_1R1C_2D_STO-4D   TYPE_of_MODEL = 8
     . MODEL=DIFFUSION_STOLLENBERG_3D   TYPE_of_MODEL = 10
     . MODEL=DIFFUSION_STOLLENBERG_4D   TYPE_of_MODEL = 15
     . MODEL=DIFFUSION_BD_2D            TYPE_of_MODEL = 13
     . MODEL=DIFFUSION_BD_3D            TYPE_of_MODEL = 14
     . MODEL=DIFFUSION_HII_2D           TYPE_of_MODEL = 9
     . MODEL=DIFFUSION_HII_1D           TYPE_of_MODEL = 12
     . MODEL=DIFFUSION_HII_nD           TYPE_of_MODEL = 16
     . MODEL=DIFFUSION_AZTECA_4D        TYPE_of_MODEL = 17
     . MODEL=DIFFUSION_AZTECA_4D_0      TYPE_of_MODEL = 18
     . MODEL=DIFFUSION_AZTECA_4D_1      TYPE_of_MODEL = 19
     . MODEL=DIFFUSION_ECOEVO_PLANTS    TYPE_of_MODEL = 20
     . MODEL=DIFFUSION_ECO_PLASMIDS     TYPE_of_MODEL = 21
     Therefore, I will make sure these are the models at work
     when the program comes to this point. 
  */
  if(P->TYPE_of_MODEL == 21 || P->TYPE_of_MODEL == 17 || P->TYPE_of_MODEL == 18 || P->TYPE_of_MODEL == 19 || P->TYPE_of_MODEL == 20 || P->TYPE_of_MODEL == 2 || P->TYPE_of_MODEL == 8 || P->TYPE_of_MODEL == 10 || P->TYPE_of_MODEL == 15 || P->TYPE_of_MODEL == 12 || P->TYPE_of_MODEL == 13 || P->TYPE_of_MODEL == 14 || P->TYPE_of_MODEL == 9 || P->TYPE_of_MODEL == 16){
    Event_Delta_Matrix_Initialization(PATCH, P);
    Event_Adjacence_List_Initialization(PATCH, P);
  }
  else{
    Print_Press_Key(1,1,"Stochastic optimization has not been implemented for this model.\n");
  }
#endif
}

void Immigration_Preassure_on_Focal_Patch_Initialization( Community ** PATCH,
							  Parameter_Model * P )
{
  /* This preassure depends on population structure across local populations */

  /* This function is required when the initial condition is set up. It calculates the
     immigration preassures on each patch caused by individuals of every disease
     status located in the neighborhood of the focal patch                   */

  /* However, this function should be only used in conjunction with:

     Immigration_Preassure_on_Focal_Patch_Update();
     Temporal_Dynamics_Update();

     which are all functions intended to optimize the algorithm. The idea is that all
     temporal rates for the different events, including immigration preassures, are only
     calculated once, at the initial condition, and then they are updated (by summing or
     substracting) according to the flow of events that occur in the system

     These two functions have not been implemented yet.
  */

  int Sp, i, j, k, n;
  double Imm_Rate;

  /* Sp is the number of variables required to define the state of a single patch */
  Sp = P->LOCAL_STATE_VARIABLES;

  for(j=0; j < Sp ; j++) {

    for(i=0; i<P->No_of_CELLS; i++) {

      Imm_Rate = 0.0;
      for(n=0; n < PATCH[i]->No_NEI; n++){
	      Imm_Rate += PATCH[i]->In_Migration_Vector[j][n]*(double)PATCH[i]->NEI[n]->n[j];
	      PATCH[i]->Imm_Rates_Preassure[j][n] = PATCH[i]->In_Migration_Vector[j][n]*(double)PATCH[i]->NEI[n]->n[j];
      }

      PATCH[i]->Total_Imm_Rate_Preassure[j] = Imm_Rate;
    }
  }
}

void Network_Structure_Inititialization (Community ** PATCH,
					                               int No_of_NEIGHBORS,
					                               int TYPE_of_NETWORK)
{
  int a, i,j,n, no, N_X, N_Y;
  int Sp;
  int i_x, j_y;
  double STEP_X, STEP_Y;
  double Total_Per_Capita_Out_Migration_Rate;

  switch ( TYPE_of_NETWORK )
  {
    case 0: /* Fully Connected Graph */

      no        = PATCH[0]->No_of_CELLS;
      Sp        = PATCH[0]->LOCAL_STATE_VARIABLES;


      for(i=0; i<no; i++){

	      PATCH[i]->center.x = gsl_rng_uniform(r) * PATCH[i]->X_DIMENSION;
	      PATCH[i]->center.y = gsl_rng_uniform(r) * PATCH[i]->Y_DIMENSION;

	        for( a=0; a<Sp; a++ ) {

	          n=0;
	          Total_Per_Capita_Out_Migration_Rate = 0.0;
	          for(j=0; j<no; j++){
	            if( i != j) {
	              PATCH[i]->NEI[n] = PATCH[j];

	              PATCH[i]->Out_Migration_Vector[a][n] = PATCH[i]->Metapop_Connectivity_Matrix[a][j][i];
	              PATCH[i]->In_Migration_Vector[a][n]  = PATCH[i]->Metapop_Connectivity_Matrix[a][i][j];

	              PATCH[i]->Patch_Connections[n] = j;

	              Total_Per_Capita_Out_Migration_Rate += PATCH[i]->Out_Migration_Vector[a][n];
	              n++;
	            }
	          }
	          PATCH[i]->Total_Per_Capita_Out_Migration_Rate[a] = Total_Per_Capita_Out_Migration_Rate;
      	  }
        PATCH[i]->No_NEI = no-1; /* All patches are connected to i */
	      assert(no-1 == n);
      }
      break;

    case 1: /* Squared Grid with Von Neuman neighborhood */

      no        = PATCH[0]->No_of_CELLS;
      Sp        = PATCH[0]->LOCAL_STATE_VARIABLES;
      N_X       = PATCH[0]->No_of_CELLS_X;
      N_Y       = PATCH[0]->No_of_CELLS_Y;

      STEP_X    = PATCH[0]->X_DIMENSION/(double)PATCH[0]->No_of_CELLS_X;
      STEP_Y    = PATCH[0]->Y_DIMENSION/(double)PATCH[0]->No_of_CELLS_Y;

      for(i=0; i<no; i++){

	      i_x = i/PATCH[i]->No_of_CELLS_X;
	      j_y = i%PATCH[i]->No_of_CELLS_X;

	      PATCH[i]->center.x = (double)j_y + 0.5*STEP_X;
	      PATCH[i]->center.y = (double)i_x + 0.5*STEP_Y;

	      Set_Von_Neumann_1st_Neighbors(PATCH, no, N_X, N_Y, i);

	      for( a=0; a<Sp; a++ ) {

	        Total_Per_Capita_Out_Migration_Rate = 0.0;
	        for(j=0; j<No_of_NEIGHBORS; j++){

	          PATCH[i]->Out_Migration_Vector[a][j] = PATCH[i]->Metapop_Connectivity_Matrix[a][i][j];
	          PATCH[i]->In_Migration_Vector[a][j]  = PATCH[i]->Metapop_Connectivity_Matrix[a][i][j];

	          Total_Per_Capita_Out_Migration_Rate += PATCH[i]->Out_Migration_Vector[a][j];
	        }

	        PATCH[i]->Total_Per_Capita_Out_Migration_Rate[a] = Total_Per_Capita_Out_Migration_Rate;
	      }
	      PATCH[i]->No_NEI = No_of_NEIGHBORS;
      }
    break;

    default:
      printf("Type of Network not yet defined!!!\n");
      printf("Allowed Network Codes are: 0 and 1\n");
      printf("TYPE of NETWORK = %d\n", TYPE_of_NETWORK);
      exit(0);
  }
}

void  Set_Von_Neumann_1st_Neighbors(Community ** PATCH, int no, int N_X, int N_Y, int i)
{
  int i_x, i_y;
  int n_x, n_y;
  int nei;

  i_x = i/N_X;
  i_y = i%N_X;

  /* Upper Neighbor */
  n_x = i_x;
  n_y = (i_y+1)%N_Y;                    /* Periodic Boundary Condition */
  nei = n_x * N_X  + n_y;
  PATCH[i]->NEI[0] = PATCH[nei];
  PATCH[i]->Patch_Connections[0] = nei;

  /* Right Neighbor */
  n_x = (i_x+1)%N_X;                    /* Periodic Boundary Condition */
  n_y = i_y;
  nei = n_x * N_X  + n_y;
  PATCH[i]->NEI[1] = PATCH[nei];
  PATCH[i]->Patch_Connections[1] = nei;

  /* Lower Neighbor */
  n_x = i_x;
  n_y = (i_y == 0) ? (N_Y-1) : (i_y-1); /* Periodic Boundary Condition */
  nei = n_x * N_X  + n_y;
  PATCH[i]->NEI[2] = PATCH[nei];
  PATCH[i]->Patch_Connections[2] = nei;

  /* Left Neighbor */
  n_x = (i_x == 0) ? (N_X-1) : (i_x-1); /* Periodic Boundary Condition */
  n_y = i_y;
  nei = n_x * N_X  + n_y;
  PATCH[i]->NEI[3] = PATCH[nei];
  PATCH[i]->Patch_Connections[3] = nei;
}

void Writing_Adjacency_List(Community ** PATCH)
{
  int i,j, a, no, n_Sp;

  no    = PATCH[0]->No_of_CELLS;
  n_Sp  = PATCH[0]->LOCAL_STATE_VARIABLES;

  for(a=0; a<n_Sp; a++){
    printf("%s %d\n", "Species", a);
    for(i=0; i<no; i++){
      printf("%s %d %s", "Local Population No", i, "is conntected to [  ");
      for(j=0; j<PATCH[i]->No_NEI; j++) {
        printf("%d ", PATCH[i]->Patch_Connections[j]);
      }
      printf("%s\n", " ].");
    }
    printf("\n");
  }

  // Print_Press_Key(1,0,".");
}

void Writing_Adjacency_List_VonNeumann(Community ** PATCH)
{
  int i,j, a, no, n_Sp;
  int i_x, i_y, N_X;
  int j_x, j_y;

  no    = PATCH[0]->No_of_CELLS;
  N_X   = PATCH[0]->No_of_CELLS_X;
  n_Sp  = PATCH[0]->LOCAL_STATE_VARIABLES;

  for(a=0; a<n_Sp; a++){
    printf("%s %d\n", "Species", a);
    for(i=0; i<no; i++){
      i_x = i/N_X;
      i_y = i%N_X;
      printf("%s %d %s (%d, %d) %s", "Local Population No", i, "located at node", i_x, i_y, "is conntected to [");
      for(j=0; j<PATCH[i]->No_NEI; j++) {
	      j_x = PATCH[i]->Patch_Connections[j]/N_X;
	      j_y = PATCH[i]->Patch_Connections[j]%N_X;
        printf("  [patch %d located at (%d, %d)]  ",
	       PATCH[i]->Patch_Connections[j], j_x, j_y);
      }
      printf("%s\n", " ].");
    }
    printf("\n");
  }

  // Print_Press_Key(1,0,".");
}

void Print_Meta_Community_Patch_System (Parameter_Table * Table)
{
  int k, Sp, Patch;

  Community ** Village = Table->Patch_System;

  Sp = Table->LOCAL_STATE_VARIABLES; /* 'The number of state variables that fully define
					the configuration of any given patch
				     */
  printf(" Total population on local populations (checking local pointers, Y and J vectors):\n");

  for(Patch=0; Patch<Table->No_of_CELLS; Patch++) {
    printf(" Patch[%d]:\t", Patch);

    for(k = 0; k < Sp; k++) {

      printf(" %s = %d ", Table->Model_Variable_Name[k + Patch*Sp],
	     Village[Patch]->n[k]);
      printf(" %s = %g ", Table->Model_Variable_Name[k + Patch*Sp],
	     Table->Vector_Model_Variables[k+Patch*Sp]);
      printf(" %s = %d ", Table->Model_Variable_Name[k + Patch*Sp],
	     Table->Vector_Model_Int_Variables[k+Patch*Sp]);
      printf("\t");
    }

    printf("\n");

  }
}

void Community_Binary_Tree_Allocation (Parameter_Table * Table, int No_of_CELLS)
{
  int i, k, No_of_LEAVES, No_of_TREE_LEVELS, No;

  /* Determine the value of No_of_LEAVES and No_of_TREE_LEVELS given that, at least, 
     there should be enough for a No_of_CELLS/GRAND_No_of_EVENTS of true active leaves 
  */ 
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

  Table->No_of_LEAVES      = No_of_LEAVES;
  Table->No_of_TREE_LEVELS = No_of_TREE_LEVELS; 

  Table->Treeroot = Binary_Tree_Allocation ( No_of_CELLS, 
                                             &(Table->Leaves), 
                                             &(Table->Parent) );
}

void Community_Binary_Tree_Initialization (Parameter_Table * Table)
{
  /* 
     This function sets up the partial sums of a previously allocated
     binary tree that will maintain the discrete probability distribution 
     always ready to be sampled. 
  */

  treenode *** Parent = Table->Parent; /* Set of Parent nodes at each level       */
  treenode ** Leaves  = Table->Leaves; /* from level 0 (root) to level n (leaves) */

  int n = Table->No_of_TREE_LEVELS; 

  Table->Treeroot = sumBinaryTree_DiscreteDistribution(Parent, Leaves, n); 
}

void Community_Priority_Queue_Tree_Allocation ( Parameter_Table * Table, 
                                                int TOTAL_GRAND_No_of_EVENTS )
{
  int l, k, No_of_LEAVES, No_of_TREE_LEVELS, No;
  double l_x; 

  /* Determine the value of No_of_LEAVES and No_of_TREE_LEVELS given that, at least, 
     there should be enough to acommodate a TOTAL_GRAND_No_of_EVENTS of true active nodes 
     over the whole tree (counting from root node, and both interval nodes and leaves). 
  */ 
  No_of_TREE_LEVELS = 0;     /* Only the root!!! */
  No_of_LEAVES      = 1;     /* The root!!!      */
 
  if ( TOTAL_GRAND_No_of_EVENTS > 1) {
    No_of_TREE_LEVELS = Calculating_No_of_TREE_LEVELS(TOTAL_GRAND_No_of_EVENTS);
    No_of_LEAVES      = power_int(2, No_of_TREE_LEVELS);
    /* No_of_TREE_LEVELS corresponds to the No of (internal) TREE LEVELS (without
       counting the leaves) and, of course, also to the label of the tree level 
       corresponding to the leaves (since the root tree level is labeled as zero).  
    */
  }

  assert( (power_int(2, No_of_TREE_LEVELS+1) - 1) > TOTAL_GRAND_No_of_EVENTS );

  Table->No_of_LEAVES      = No_of_LEAVES;
  Table->No_of_TREE_LEVELS = No_of_TREE_LEVELS; /* No of (internal) TREE LEVELS (without counting 
                                                   the final outer LEAVE level)
                                                */
  Table->Tree_Node_Index = malloc(TOTAL_GRAND_No_of_EVENTS * sizeof(treenode *));
                                           
  Table->Treeroot = Binary_Tree_Allocation ( No_of_LEAVES, 
                                             &(Table->Leaves), 
                                             &(Table->Parent) );
}

void Community_Priority_Queue_Tree_Initialization (Parameter_Table * Table)
{
  /* 
     This function sets up the priority queu that will maintain the minimum
     time value at the root level (in spite of changes in any internal 
     tree node). See treenode.c library (bubbling algorithm).
  */
  int i, N, n, n_0, m, n_1; 
  double Next_Time; 

  treenode *** Parent = Table->Parent; /* Set of Parent nodes at each level       */
  treenode ** Leaves  = Table->Leaves; /* from level 0 (root) to level n (leaves) */
  treenode ** Priority = Table->Tree_Node_Index; 

  double   * Vector   = Table->T->Vector_of_Rates;

  for(i=0; i<Table->TOTAL_GRAND_No_of_EVENTS; i++) {
    if(Vector[i] > 0.0)
      Next_Time = -1/Vector[i] * log (RANDOM);
    else
      Next_Time = INFINITY; 

    Priority_Queu_Insert_Value(i, Next_Time, Table->No_of_TREE_LEVELS, 
                               Priority, Parent, Leaves); 
  } 
  /* Those leaves that are not in use should be made: 
        . index >= Table->TOTAL_No_of_EVENTS 
        . value  = INFINITY 
     because they should represent "impossible events" that will never 
     bubble up to the root!!! 
  */
  /* Calculationg Leaves not in use */
  N = Table->TOTAL_GRAND_No_of_EVENTS;
  n = Table->No_of_TREE_LEVELS;
  n_0 = N - power_int(2, n) + 1;
  n_1 = power_int(2, n); 
  m   = N; 
  for(i=n_0; i<n_1; i++) {
    Table->Leaves[i]->index = m++;  
    Table->Leaves[i]->value = INFINITY; 
  } 
}

#if defined DIFFUSION_ECO_PLASMIDS
void Community_Plasmids_Initialization (Community ** PATCH, Parameter_Model * P)
{
  int i, j;
 
  Parameter_Table * Table = (Parameter_Table *)P->Table; 

  for (j = 0; j < Table->No_of_CELLS; j++) {

    for(i=0; i < Table->No_of_PLASMIDS; i++){

      PATCH[j]->Local_Plasmid_Population[i]->ID = i;

      PATCH[j]->Local_Plasmid_Population[i]->n  = 0; 
        
      PATCH[j]->Local_Plasmid_Population[i]->Cost = Table->Alpha_C_0;     /* Alpha_C_0 = 1.0: Full Cost: No reproduction!!! */

      PATCH[j]->Local_Plasmid_Population[i]->Resistance = Table->Nu_C_0;  /* Nu_C_0 = 1.0: Full Resistance */

      PATCH[j]->Local_Plasmid_Population[i]->Compatibility_Profile = Table->Plasmid_Compatibility_Indeces[i];
    }
  }
}

void Community_Strains_Initialization (Community ** PATCH, Parameter_Model * P )
{
  int i, j, i_Strain, k_Profile;
  int Sp;
  double COST, RESISTANCE; 

  Sp    = P->LOCAL_STATE_VARIABLES;  /* Total Number of Species 
					                              potentially coexisting locally,
					                              and, therefore, also the
					                              Total Number of State Variables
					                              fully determining the local state
  				                            */ 
  Parameter_Table * Table = (Parameter_Table *)P->Table; 

  for (j = 0; j < Table->No_of_CELLS; j++) {

    for(i=0; i < Sp; i++){

      Calculate_Strain_and_Profile(Table, i, &i_Strain, &k_Profile);

      PATCH[j]->Local_Strain_Population[i]->ID = i;

      PATCH[j]->Local_Strain_Population[i]->n  = 0;  

      PATCH[j]->Local_Strain_Population[i]->Death_0 = Table->Death_AP[i];

      PATCH[j]->Local_Strain_Population[i]->Death_1 = Table->Death_R_1; 

      PATCH[j]->Local_Strain_Population[i]->Beta = Table->Beta_R[i];

      PATCH[j]->Local_Strain_Population[i]->Segregation_Error = Table->Segregation_Error[i];

      PATCH[j]->Local_Strain_Population[i]->Gamma = Table->Eta_RP[i];

      PATCH[j]->Local_Strain_Population[i]->Mu_0 = Table->Mu_RP[i];

      PATCH[j]->Local_Strain_Population[i]->Competition_List = Table->Competition_List_Indeces[i]; 

      PATCH[j]->Local_Strain_Population[i]->Conjugation_List = Table->Conjugation_List_Indeces[i];

      PATCH[j]->Local_Strain_Population[i]->Recipient_List = Table->Recipient_List_Indeces[i];

      PATCH[j]->Local_Strain_Population[i]->Donor_List = Table->Donor_List_Indeces[i];

      PATCH[j]->Local_Strain_Population[i]->Profile = Table->Strain_Profiles[i_Strain][k_Profile];
    }
  }
}
#endif