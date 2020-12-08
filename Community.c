#include <MODEL.h>

/* This functions allocate, initialize and free a number of local communities,
   which make up our total patch system or metapopulation */

void Community_Allocation ( Community ** PATCH, Parameter_Model * P )
{
  int i, j, a, K;
  int no, Sp;

  Sp    = P->No_of_SPECIES; 
  no    = P->No_of_CELLS;

  for(i=0; i<no; i++){
    PATCH[i] = (Community *)calloc( 1, sizeof(Community) );

    PATCH[i]->n = (int *)calloc(Sp, sizeof( int ));

    PATCH[i]->rate = (double *)calloc(P->TOTAL_No_of_EVENTS, sizeof( double ));
    PATCH[i]->rToI = (double *)calloc(P->TOTAL_No_of_EVENTS, sizeof( double ));

    PATCH[i]->Patch_Connections = (int *)calloc(P->No_of_NEIGHBORS, sizeof( double ));

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
      PATCH[i]->Imm_Rates_Preassure[a] = (double *)calloc(P->No_of_NEIGHBORS,
							   sizeof( double ));
  }
}

void Community_Free (Community ** PATCH, Parameter_Model * P)
{
  int Sp, K, i, j, a;

  Sp  = P->No_of_SPECIES; /* Ex: 11 times 4 */

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

    free(PATCH[i]);
  }

  free( PATCH );
  /*   END: Patch Total Destruction */
}

void Community_Initialization (Community ** PATCH,
			       Parameter_Model * P )
{
  int i, j, Sp, no;
  int i_x, j_y; 
  double STEP;
  
  Sp  = P->No_of_SPECIES;
  no  = P->No_of_CELLS;

  STEP = 1.0; 
    
  for(i=0; i<no; i++){

    i_x = i/P->No_of_CELLS_X;
    j_y = i%P->No_of_CELLS_X;
    
    PATCH[i]->No_of_SPECIES = P->No_of_SPECIES;

    PATCH[i]->No_of_CELLS   = P->No_of_CELLS;
    PATCH[i]->No_of_CELLS_X = P->No_of_CELLS_X;
    PATCH[i]->No_of_CELLS_Y = P->No_of_CELLS_Y;
    PATCH[i]->X_DIMENSION = (double)P->No_of_CELLS_Y;
    PATCH[i]->Y_DIMENSION = (double)P->No_of_CELLS_X;
    
    
    PATCH[i]->center.x = (double)j_y + 0.5*STEP;
    PATCH[i]->center.y = (double)i_x + 0.5*STEP;

    PATCH[i]->no_VARIABLES  = Sp;
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

  Immigration_Preassure_on_Focal_Patch_Initialization( PATCH, P );

  Writing_Adjacency_List(PATCH);

  if (P->TYPE_of_NETWORK == 1) Writing_Adjacency_List_VonNeumann(PATCH);
                               
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
  Sp = P->No_of_SPECIES;

  for(i=0; i<P->No_of_CELLS; i++) {
    for(j=0; j < Sp ; j++) {

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

  double Total_Per_Capita_Out_Migration_Rate = 0.0;

  switch ( TYPE_of_NETWORK )
    {

    case 0: /* Fully Connected Graph */

      no        = PATCH[0]->No_of_CELLS;
      Sp        = PATCH[0]->No_of_SPECIES;

      for( a=0; a<Sp; a++ ) {
       
	for(i=0; i<no; i++){
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
	  PATCH[i]->No_NEI = n; /* All patches are connected to i */
	  PATCH[i]->Total_Per_Capita_Out_Migration_Rate[a] = Total_Per_Capita_Out_Migration_Rate;
	}
      }
      break;

    case 1: /* Squared Grid with Von Neuman neighborhood */

      no        = PATCH[0]->No_of_CELLS;
      Sp        = PATCH[0]->No_of_SPECIES;
      N_X       = PATCH[0]->No_of_CELLS_X;
      N_Y       = PATCH[0]->No_of_CELLS_Y;

      for( a=0; a<Sp; a++ ) {	
	for(i=0; i<no; i++){
	  
	  Set_Von_Neumann_1st_Neighbors(PATCH, no, N_X, N_Y, i);

	  Total_Per_Capita_Out_Migration_Rate = 0.0;
	  for(j=0; j<No_of_NEIGHBORS; j++){
	    PATCH[i]->Out_Migration_Vector[a][j] = PATCH[i]->Metapop_Connectivity_Matrix[a][i][j];
	    PATCH[i]->In_Migration_Vector[a][j]  = PATCH[i]->Metapop_Connectivity_Matrix[a][i][j];

	    Total_Per_Capita_Out_Migration_Rate += PATCH[i]->Out_Migration_Vector[a][j];
	  }

	  PATCH[i]->No_NEI = No_of_NEIGHBORS;
	  PATCH[i]->Total_Per_Capita_Out_Migration_Rate[a] = Total_Per_Capita_Out_Migration_Rate;
	}
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
  n_Sp  = PATCH[0]->No_of_SPECIES;

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

  // Press_Key();
}

void Writing_Adjacency_List_VonNeumann(Community ** PATCH)
{
  int i,j, a, no, n_Sp;
  int i_x, i_y, N_X;
  int j_x, j_y;

  no    = PATCH[0]->No_of_CELLS;
  N_X   = PATCH[0]->No_of_CELLS_X;
  n_Sp  = PATCH[0]->No_of_SPECIES;

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

  // Press_Key();
}
