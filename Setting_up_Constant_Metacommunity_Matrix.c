#include <MODEL.h>

void Setting_up_Constant_Metapopulation_Connectivity_Matrix (Parameter_Table * Table) 
{
	/* This function sets up the connectity matrix for a regular network (each node has 
	   the same number of neighbors, Table->No_of_NEIGHBORS). In fact, different species
	   in the system can differ in movility (different Mu rates). The number of different
	   species in the system is given by Table->LOCAL_STATE_VARIABLES. 
	   The function sets up a differetn Adjancency Matrix for every species in the system.

	   main() (main.c) 
          --->  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( ... ) (Parameter_Table.c)
            ---> Setting_up_Constant_Metapopulation_Connectivity_Matrix (...)                      
	*/
    int i, j, a;

    if( Table->TYPE_of_MODEL == 0 || Table->TYPE_of_MODEL == 1 )
        for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++) 
	      for(i=0; i<Table->No_of_CELLS; i++)
	        for(j=0; j<Table->No_of_NEIGHBORS; j++)
	          Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;
    
    else if (Table->TYPE_of_MODEL == 2 || Table->TYPE_of_MODEL == 10 || Table->TYPE_of_MODEL == 13)
        for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
	      if(a == 0)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;
	      else if (a == 1)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu_C;
          else if (a == 2)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	          Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0 * Table->Mu_C;
	      else 
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0;
    
    else if (Table->TYPE_of_MODEL == 15 || Table->TYPE_of_MODEL == 16)  
    /* DIFFUSION_STOLLENBERG_4D or DIFFUSION_HII_nD */
      for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
	      if(a == 0)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;
	      else if (a == 1)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0; 
        else if (a == 2)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu_C;
	      else 
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0;

    else if (Table->TYPE_of_MODEL == 17 || Table->TYPE_of_MODEL == 18)  
    /* DIFFUSION_AZTECA_4D or DIFFUSION_AZTECA_4D_0*/
      for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
	      if(a == 0)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;   /* Workers */
	      else if (a == 1)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0;         /* Queens  */
        else if (a == 2)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu_C; /* Flies   */
	      else 
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0;         /* Parasitized Workers */ 
    
    else if (Table->TYPE_of_MODEL == 19)  
    /* DIFFUSION_AZTECA_4D_1*/
      for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
	      if(a == 0)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;   /* Workers */
	      else if (a == 1)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Lambda_C_1; /* Queens  */
          else if (a == 2)
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu_C; /* Flies   */
	      else 
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0.0;  /* Parasitized Workers */ 

    else if (Table->TYPE_of_MODEL == 20)  
    /* DIFFUSION_ECOEVO_PLANTS */
      for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
	      if(a%2 == 0)                                           /* (a%2 == 0, this is, resource propagules) */
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu;   /* Propagues, RP, move at rate Mu      */
	      else                                                   /* (a%2 == 1, this is, established adult plants) */
	        for(i=0; i<Table->No_of_CELLS; i++)
	          for(j=0; j<Table->No_of_NEIGHBORS; j++)
	            Table->Metapop_Connectivity_Matrix[a][i][j] = 0;           /* Adult Plants, AP, on't move themseleves */
    
    else if (Table->TYPE_of_MODEL == 21)  
    /* DIFFUSION_ECOEVO_PLASMIDS */
      for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
	      for(i=0; i<Table->No_of_CELLS; i++)
	        for(j=0; j<Table->No_of_NEIGHBORS; j++)
	          Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu_RP[a];  /* All bacterial strains move at 
                                                                               the same rate Mu */
                                                                              /* This can be reset according to a 
                                                                                 strain diffusion vector 
                                                                              */
    else if (Table->TYPE_of_MODEL == 22)  
   /* DIFFUSION_ECO_1B1P */
      for(a=0; a<Table->LOCAL_STATE_VARIABLES; a++)
        for(i=0; i<Table->No_of_CELLS; i++)
          for(j=0; j<Table->No_of_NEIGHBORS; j++)
            Table->Metapop_Connectivity_Matrix[a][i][j] = Table->Mu_RP[a];   /* The plasmid free bacterial strain 
                                                                                moves at rate Mu  */
                                                                             /* The plasmid bearing bacterial strain 
                                                                                moves at rate Mu_C */
    else{
        fprintf(stderr, " TYPE of MODEL (%d) not defined (at Parameter_Table.c)\n",
                Table->TYPE_of_MODEL);
        printf(" This model is not ready for multi-patch dynamics\n");
        assert(Table->No_of_CELLS == 1); 
        exit(0); 
    }
}                                                                         