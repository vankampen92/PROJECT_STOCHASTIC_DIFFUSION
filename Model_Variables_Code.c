#include <MODEL.h>

void Model_Variables_Code_into_Parameter_Table (Parameter_Table * Table)
{
  int i, j, n;
  int TYPE_of_MODEL;

  TYPE_of_MODEL = Table->TYPE_of_MODEL;

  switch( TYPE_of_MODEL )
    {

    case 0: /* DIFFUSION * * * * * * * * * * * * * * * * * * * * * * */

      /* Number of events that can occur to a single Species: */
      Table->No_of_EVENTS       = 1;  /* (Only Diffusion)         */
      Table->TOTAL_No_of_EVENTS = Table->No_of_EVENTS * Table->No_of_RESOURCES;
      Table->LOCAL_STATE_VARIABLES = Table->No_of_RESOURCES;
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */

      n = 0; 
      Table->Index[n++] = 0; /* Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */  

      Table->TOTAL_No_of_MODEL_PARAMETERS = n; 
      
      break;

    case 1: /* DIFFUSION_S_RESOURCES * * * * * * * * * * * * * * * * * * * * * * */
      
      /* Number of events that can occur to a single Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = Table->No_of_EVENTS * Table->No_of_RESOURCES;
      Table->LOCAL_STATE_VARIABLES = Table->No_of_RESOURCES;
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */

      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      if ( Table->No_of_RESOURCES > 0 ) {
	Table->Index[n++] = 6;
	Table->Index[n++] = 7; 
      }
      if ( Table->No_of_RESOURCES > 1 ) {
	Table->Index[n++] = 8;
	Table->Index[n++] = 9; 
      }
      Table->Index[n++]   = 10; /* No_of_RESOURCES */ 
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      
      break;

    case 2: /* DIFFUSION_1R1C * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 5;
      Table->LOCAL_STATE_VARIABLES = 4; /* 1 R + 1 C + 1 D + 1 T        */
                                        /* D \equiv RC and T \equiv CRC */
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->A = 1; Table->RA = 2; Table->ARA = 3; 
      
      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Resource Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      Table->Index[n++] = 6; /* External Immigration Rate (0) */
      Table->Index[n++] = 7; /* Death Rate (0) */
      
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
      Table->Index[n++]   = 18; /* Trimer Formation Rate */
      Table->Index[n++]   = 19; /* Trimer Destruction Rate */
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 3: /* DIFFUSION_1RnC_E * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 7;
      Table->LOCAL_STATE_VARIABLES = 1 + 2*Table->N_E + Table->N_E*Table->N_E;
                                   /* 1R + N_E + N_E + N_E*N_E    */
    
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */

      Table->SUM_LOCAL_STATE_VARIABLES = 5; 
      Table->R = 0; Table->A = 1; Table->A_R = 2; Table->RA = 3; Table->ARA = 4;
      /* A_R: Index for the mature population contributing to reproduction */
      
      for(i=0; i<Table->N_E; i++) {
	Table->A_P[i]  = Table->R + 1 + i;
	Table->RA_P[i] = Table->R + Table->N_E + 1 + i ; 
      }
      for(i=0; i<Table->N_E; i++)
	for(j=0; j<Table->N_E; j++) 
	  Table->ARA_P[i][j] = i + j*Table->N_E + Table->R + 2*Table->N_E + 1; 

      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Resource Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      Table->Index[n++] = 6; /* External Immigration Rate (0) */
      Table->Index[n++] = 7; /* Death Rate (0) */
      
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
      Table->Index[n++]   = 18; /* Trimer Formation Rate */
      Table->Index[n++]   = 19; /* Trimer Destruction Rate */
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */

      Table->Index[n++]   = 21; /* Number of Energy Levels */
      Table->Index[n++]   = 22; /* Fecundity: Number of Offspring Individuals */
      Table->Index[n++]   = 23; /* Energy Level at Maturity  */
      Table->Index[n++]   = 24; /* Consummer Reproduction Rate */
      Table->Index[n++]   = 25; /* 2* k_E is the resourse value in energy units */
      Table->Index[n++]   = 26; /* Energy loss rate for maintenance */
      Table->Index[n++]   = 27; /* Cooperation probability 1st position in the triplet */ 
      Table->Index[n++]   = 28; /* Cooperation probability 2on position in the triplet */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break; 

    case 4: /* DIFFUSION_1R1C_2D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 5;
      Table->LOCAL_STATE_VARIABLES = 2; /* 1 R + 1 C */
                                        
      assert(Table->No_of_RESOURCES == 1);

      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
        for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
          n++;

      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->A = 1; 

      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Resource Movement Rate */
      Table->Index[n++] = 5; /* No_of_RESOURCES */

      if ( Table->No_of_RESOURCES > 0 ) {
        Table->Index[n++] = 6; /* External Immigration Rate (0) */
        Table->Index[n++] = 7; /* Death Rate (0) */
      }
      if ( Table->No_of_RESOURCES > 1 ) {
        Table->Index[n++] = 8;
        Table->Index[n++] = 9;
      }
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */

      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */

      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */

      Table->Index[n++]   = 20; /* Consumer Movement Rate  */

      Table->Index[n++]   = 26; /* Energy Loss Rate (here, fraction of 
				   conumed resourced that are NOT transformed 
				   into new consumers 
				*/

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break; 

    case 5: /* DIFFUSION_DRAG * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species:                      */
      Table->No_of_EVENTS       = 3;    /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 3 * Table->No_of_EVENTS + 5; // 5? Stochastic Dynamics Not Yet Implemented 
      Table->LOCAL_STATE_VARIABLES = 3; /* V + R + G */
                                        
      assert(Table->No_of_RESOURCES == 1);

      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
        for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
          n++;

      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->V = 0;  Table->R = 1; Table->G = 2; 

      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0;  /* Resource Movement Rate     */

      Table->Index[n++] = 6;  /* External Immigration Rate (V)    */
      Table->Index[n++] = 7;  /* Per Capita Death Rate (V)        */
      Table->Index[n++] = 10; /* Vegetation Carrying Capacity (V) */
      Table->Index[n++] = 11; /* Vegetation Local Growth Rate (V) */
      
      Table->Index[n++] = 12; /* External Immigration Rate (R) */
      Table->Index[n++] = 13; /* Rat Per Capita Death Rate (R) */  

      Table->Index[n++] = 14; /* External Immigration Rate (G) */
      Table->Index[n++] = 15; /* Rat Per Capita Death Rate (G) */  
      
      Table->Index[n++] = 16; /* Rat Exploration/Attack Rate (R) */
      
      Table->Index[n++] = 20; /* Consumer Movement Rate (G / R)  */

      Table->Index[n++] = 24; /* Gull Per Capita Reproduction Rate (G) */

      Table->Index[n++] = 26; /* Energy Loss/Efficient Rate (here, fraction of 
				 conumed resources that are NOT transformed 
				 into new consumers) */
      Table->Index[n++] = 27; /* Guano Conversion Factor (G) */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break; 

    case 6: /* DIFFUSION_VRG * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species:                      */
      Table->No_of_EVENTS       = 3;    /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 3 * Table->No_of_EVENTS + 5; // 5? Stochastic Dynamics Not Yet Implemented 
      Table->LOCAL_STATE_VARIABLES = 3; /* V + R + G */
                                        
      assert(Table->No_of_RESOURCES == 1);

      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
        for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
          n++;

      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->V = 0;  Table->R = 1; Table->G = 2; 

      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0;  /* Resource Movement Rate     */

      Table->Index[n++] = 6;  /* External Immigration Rate (V)    */
      Table->Index[n++] = 7;  /* Per Capita Death Rate (V)        */
      Table->Index[n++] = 10; /* Vegetation Carrying Capacity (V) */
      Table->Index[n++] = 11; /* Vegetation Local Growth Rate (V) */
      
      Table->Index[n++] = 12; /* External Immigration Rate (R) */
      Table->Index[n++] = 13; /* Rat Per Capita Death Rate (R) */  

      Table->Index[n++] = 14; /* External Immigration Rate (G) */
      Table->Index[n++] = 15; /* Rat Per Capita Death Rate (G) */  
      
      Table->Index[n++] = 16; /* Rat Exploration/Attack Rate (R) */
      
      Table->Index[n++] = 20; /* Consumer Movement Rate (G / R)  */

      Table->Index[n++] = 24; /* Gull Per Capita Reproduction Rate (G) */

      Table->Index[n++] = 26; /* Energy Loss/Efficient Rate (here, fraction of 
				 conumed resources that are NOT transformed 
				 into new consumers) */
      Table->Index[n++] = 27; /* Guano Conversion Factor (G) */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break; 

    case 7: /* DIFFUSION_MR * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 5;
      Table->LOCAL_STATE_VARIABLES = 3; /* 1 R + 1 C + 1 D + 1 T        */
                                        /* D \equiv RC and T \equiv CRC */
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->A = 1; Table->RA = 2; 
      
      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Resource Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      if ( Table->No_of_RESOURCES > 0 ) {
	Table->Index[n++] = 6; /* External Immigration Rate (0) */
	Table->Index[n++] = 7; /* Death Rate (0) */
      }
      if ( Table->No_of_RESOURCES > 1 ) {
	Table->Index[n++] = 8;  
	Table->Index[n++] = 9; 
      }
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 8: /* DIFFUSION_1R1C_2D_STO_4D * * * * * * * * * * * * * * * * * * * * * * */

      assert(Table->Chi_C_0 == 0.0); 
      
      /* No_of_EVENTS: Common events that can occur to every Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 5;
      Table->LOCAL_STATE_VARIABLES = 4; /* 1 R + 1 C + 1 D + 1 T        */
                                        /* D \equiv RC and T \equiv CRC */
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
      
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->A = 1; Table->RA = 2; Table->ARA = 3; 
      
      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 0; /* Resource Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 
      
      if ( Table->No_of_RESOURCES > 0 ) {
	Table->Index[n++] = 6; /* External Immigration Rate (0) */
	Table->Index[n++] = 7; /* Death Rate (0) */
      }
      if ( Table->No_of_RESOURCES > 1 ) {
	Table->Index[n++] = 8;  
	Table->Index[n++] = 9; 
      }
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */

      Table->Index[n++]   = 26; /* Energy Loss Rate (here, fraction of 
				   conumed resourced that are NOT transformed 
				   into new consumers 
				*/
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 9: /* DIFFUSION_HII_2D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species: */
      Table->No_of_EVENTS       = 0;  
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 4;
      /* Four total events: Diffusion + Immigration + Attack + Handling */
      Table->LOCAL_STATE_VARIABLES = 2; /* 1 A + 1 RA        */
                                        /* Only free consumers (A) and
					   handling consumers  (RA) 
					*/
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->A = 0; Table->RA = 1; Table->ARA = 0; 
      
      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */

      Table->Index[n++]   = 20; /* Consumer Movement Rate  */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 10: /* DIFFUSION_STOLLENBERG_3D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every Species: */
      Table->No_of_EVENTS       = 3;  /* (Only Diffusion + External Immigration + Death) 
				          R and A can undergo these three same processes 
				      */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 5;
      Table->LOCAL_STATE_VARIABLES = 3; /* 1 R + 1 C + 1 D                 */
                                        /* D \equiv RC (handling consumer) */
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class */
      Table->R = 0;  Table->A = 1; Table->RA = 2; 

      /* List of (Potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
      Table->Index[n++] = 0; /* Resource Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      Table->Index[n++] = 6; /* External Immigration Rate (0) */
      Table->Index[n++] = 7; /* Death Rate (0) */
       
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */

      Table->Index[n++]   = 24; /* Consumer Local Growth Rate */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 11: /* DIFFUSION_HII_AC_2D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species:  */
      Table->No_of_EVENTS       = 0;  /* No of common events           */
      Table->TOTAL_No_of_EVENTS = 4;  /* Diffusion + Immigration + Attack + Handling */
      Table->LOCAL_STATE_VARIABLES = 2; /* 1 A + 1 AC        */
                                        /* Only free consumers (A) and
					   Accumulated Resource Consumed  (AC)
					   y[RA] will be calculated as A_0 - y[A]
					*/
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->AC = 0; Table->A = 1; Table->RA = 0; Table->ARA = 0; 
      
      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
     
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 12: /* DIFFUSION_HII_1D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species:  */
      Table->No_of_EVENTS       = 0;  /* No of common events           */
      Table->TOTAL_No_of_EVENTS = 4;  /* Diffusion + Ext Immigration + Attack + Handling */
      Table->LOCAL_STATE_VARIABLES = 1; /* 1 A         */
                                        /* Only free consumers (A).				   
					   y[RA] will be calculated as A_0 - y[A]
					*/
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->AC = 0; Table->A = 0; Table->RA = 0; Table->ARA = 0; 
      
      /* List of (Potentially searcheable) model parameters */
      n = 0;
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
     
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 13: /* DIFFUSION_BD_2D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species:  */
      Table->No_of_EVENTS       = 0;  /* No of common events           */
      Table->TOTAL_No_of_EVENTS = 6;  /* OutDiffusion(A) + InMigration(A) + Attack + Handling  + 
				         Triplet formation + Triplet decomposition 
				      */
      Table->LOCAL_STATE_VARIABLES = 2; /* 1 A + 1 RA        */
                                        /* Only free consumers (A) and
					   handling consumers (RA).
					   y[ARA] will be calculated as 
					   0.5*(A_0 - y[A] -y[RA])
					*/
      
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->AC = 0; Table->A = 0; Table->RA = 1; Table->ARA = 0; 
      
      /* List of (Potentially searcheable) model parameters */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */

      Table->Index[n++]   = 18; /* Trimer Formation Rate (Chi_C_0)   */
      Table->Index[n++]   = 19; /* Trimer Destruction Rate (Eta_C_0) */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  (A) */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 14: /* DIFFUSION_BD_3D * * * * * * * * * * * * * * * * * * * * * * */

      /* No_of_EVENTS: Common events that can occur to every Species:  */
      Table->No_of_EVENTS       = 0;  /* No of common events           */
      Table->TOTAL_No_of_EVENTS = 6;  /* OutDiffusion(A) + InMigration(A) + Attack + Handling  + 
				         Triplet formation + Triplet decomposition 
				      */
      Table->LOCAL_STATE_VARIABLES = 3; /* 1 A + 1 RA + 1 ARA        */
                                        
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */
      Table->R = 0;  Table->AC = 0; Table->A = 0; Table->RA = 1; Table->ARA = 2; 
      
      /* List of (Potentially searcheable) model parameters */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */

      Table->Index[n++]   = 18; /* Trimer Formation Rate (Chi_C_0)   */
      Table->Index[n++]   = 19; /* Trimer Destruction Rate (Eta_C_0) */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  (A) */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;
      
    default:
      printf(" This TYPE_of_MODEL (%d) code is not defined.\n", TYPE_of_MODEL);
      printf(" Models (0 to 10): Check input argument list!!!\n");
      exit(0);
   }
  /* Conventionally, the last label in the argument list of

     Model_Variables_Code_into_Parameter_Table (...),

     (*K), should be the label of the last model state variable.
     Then ( * K) + 1 becomes de total number of dynamic variables.
  */
  printf("Total No of MODEL PARAMETERS : %d\n", Table->TOTAL_No_of_MODEL_PARAMETERS);
  Press_Key(); 
}

void Model_Variables_Code_into_Parameter_Model (Parameter_Model * P)
{
  int i, j, n;
  int TYPE_of_MODEL;

  TYPE_of_MODEL = P->TYPE_of_MODEL;

  switch( TYPE_of_MODEL )
    {

    case 0: /* DIFFUSION * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 1: /* DIFFUSION_S_RESOURCES * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      P->K   = n-1;     /* Label last class            */

      break;

    case 2: /* DIFFUSION_1R1C * * * * * * * * * * * * * * * * * * * * * * */
      
      assert(P->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
	    
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 3: /* DIFFUSION_1RnC_E * * * * * * * * * * * * * * * * * * * * * * */
    
      assert(P->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;	    
      /* Conventions */
      P->K   = n-1;     /* Label last class            */

      break; 

    case 4: /* DIFFUSION_1R1C_2D * * * * * * * * * * * * * * * * * * * * * * */
                                  
      assert(P->No_of_RESOURCES == 1);

      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
        for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
          n++;

      /* Conventions */
      P->K   = n-1;     /* Label last class            */

      break;

    case 5: /* DIFFUSION_DRAG * * * * * * * * * * * * * * * * * * * * * * */
                                  
      assert(P->No_of_RESOURCES == 1);

      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
        for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
          n++;

      /* Conventions */
      P->K   = n-1;     /* Label last class            */

      break; 

    case 6: /* DIFFUSION_VRG * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 7: /* DIFFUSION_MR* * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 8: /* DIFFUSION_1R1C_2D_STO_4D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 9: /* DIFFUSION_HII_2D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 10: /* DIFFUSION_STOLLENBERG_3D * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 11: /* DIFFUSION_HII_AC_2D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 12: /* DIFFUSION_HII_1D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 13: /* DIFFUSION_BD_2D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break; 

    case 14: /* DIFFUSION_BD_3D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	  n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;


    default:
      printf(" This TYPE_of_MODEL (%d) code is not defined.\n", TYPE_of_MODEL);
      printf(" Check input argument list\n");
      exit(0);
    }
  
  /* Conventionally, the last label in the argument list of
     (*K), should be the label of the last model state variable.
     Then ( * K) + 1 becomes de total number of dynamic variables.
  */
}
