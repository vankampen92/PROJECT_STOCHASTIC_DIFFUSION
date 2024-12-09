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

    case 15: /* DIFFUSION_STOLLENBERG_4D * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every Species: */
      Table->No_of_EVENTS = 3;  /* (Only Diffusion + External Immigration + Death) 
				                            P and A can undergo these three same processes 
				                        */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 7;
      Table->LOCAL_STATE_VARIABLES = 4; /* 1 RP + 1 R + 1 C + 1 D          */
                                        /* D \equiv RA (handling consumer) */
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class */
      Table->RP = 0; Table->R = 1;  Table->A = 2; Table->RA = 3; 

      /* List of the 14 (potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
      Table->Index[n++] = 0; /* Propagule Movement Rate */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      Table->Index[n++] = 6; /* External Immigration Rate (0) */
      Table->Index[n++] = 7; /* Death Rate (0) */
      Table->Index[n++] = 9; /* Death Rate (1) */
       
      Table->Index[n++]   = 10; /* Resource Carrying Capacity */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */

      Table->Index[n++]   = 12; /* Consumer External Immigration Rate */
      Table->Index[n++]   = 13; /* Consumer Death Rate */
      
      Table->Index[n++]   = 16; /* Consumer Attack Rate */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */

      Table->Index[n++]   = 24; /* Consumer Local Growth Rate */

      Table->Index[n++]   = 29; /* Propagule Establishment Rate */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 16: /* DIFFUSION_HII_nD * * * * * * * * * * * * * * * * * * * * * * */
	    /* Number of events that can occur to a single Species: */
      Table->No_of_EVENTS       = 2;  /* (Attack + Handling) */
      Table->TOTAL_No_of_EVENTS = Table->No_of_EVENTS * Table->No_of_RESOURCES + 2;
      /* The last two events correspond to diffusion or jump of a free consumers 
         from the focal cell to a neighboring cell and external immigration of 
         free consumers respectively.  
      */
      Table->LOCAL_STATE_VARIABLES = Table->No_of_RESOURCES;
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
        /* The local state variables describe the number of handling consumers for each 
	      resource type */
      
      /* Conventions */
      Table->K   = n-1;     /* Label last class            */	   
      Table->RP = 0; Table->R = 0;  Table->A = 0; Table->RA = 0; 

      /* List of (Potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
    
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 
       
      Table->Index[n++] = 10; /* Resource Carrying Capacity */ 
      
      Table->Index[n++] = 16; /* Consumer Attack Rate (on average), Alpha_C_0*/

      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time (on average) on Resource 0 */

      Table->Index[n++]   = 6;  /* Nu = 1/Tau, Tau, Handling Time (on average) on Resource 1 */

      Table->Index[n++]   = 8;  /* Nu = 1/Tau, Tau, Handling Time (on average) on Resource 2 */
      
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 17: /* AZTECA_4D * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every/some of the species: */
      Table->No_of_EVENTS = 3;  /* (Only Diffusion + External Immigration + Death) 
				                            W and F can undergo these three same processes 
				                        */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 6;
      Table->LOCAL_STATE_VARIABLES = 4; /* 1 W + 1 Q + 1 F + 1 WF           */
                                        /* WF \equiv RA (handling consumer or
                                           fly larva developing or )     
                                        */
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class */
      Table->W = 0; Table->Q = 1;  Table->F = 2; Table->WF = 3; 

      /* List of the 14 (potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
      Table->Index[n++] = 0; /* Worker Movement Rate */ /* mu_R */ 
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 

      Table->Index[n++] = 6; /* External Immigration Rate (0) */ /* W */
      Table->Index[n++] = 7; /* Death Rate (0) */ /* W */
      Table->Index[n++] = 9; /* Death Rate (1) */ /* Q */
       
      Table->Index[n++]   = 10; /* Resource Carrying Capacity per Nest */ /* k */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */ /* Beta_R */ /* Q */

      Table->Index[n++]   = 12; /* External Immigration Rate */ /* F */ /* Lambda_C_0 */
      Table->Index[n++]   = 13; /* Consumer Death Rate */ /* F */ /* Delta_C_0 */

      Table->Index[n++]   = 14; /* Max Number of Nests per Patch */ /* Lambda_C_1 */    
      
      Table->Index[n++]   = 15; /* WF Death Rate */ /* Delta_C_1 *//* Delta_WF */

      Table->Index[n++]   = 16; /* Consumer Attack Rate */ /* Alpha */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */ /* Nu */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */ /* mu_C */ /* mu_F */

      Table->Index[n++]   = 29; /* Establishment Rate  (nest starting) *//* W ---> Q */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 18: /* AZTECA_4D_0 * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every/some of the species: */
      Table->No_of_EVENTS = 3;  /* (Only Diffusion + External Immigration + Death) 
				                            W and F can undergo these three same processes 
				                        */
      Table->TOTAL_No_of_EVENTS = 2 * Table->No_of_EVENTS + 6;
      Table->LOCAL_STATE_VARIABLES = 4; /* 1 W + 1 Q + 1 F + 1 WF           */
                                        /* WF \equiv RA (handling consumer or
                                           fly larva developing or )     
                                        */
      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class */
      Table->W = 0; Table->Q = 1;  Table->F = 2; Table->WF = 3; 

      /* List of the 14 (potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
  
      Table->Index[n++] = 5; /* No_of_RESOURCES */ 
      
      Table->Index[n++] = 0; /* Worker Movement Rate */ /* mu_R */ 
      
      Table->Index[n++] = 6; /* External Immigration Rate (0) */ /* W */
      
      Table->Index[n++] = 7; /* Death Rate (0) */ /* W */   /* Delta_R_0 */
      Table->Index[n++] = 9; /* Death Rate (1) */ /* Q */   /* Delta_R_1 */
       
      Table->Index[n++]   = 10; /* Next Carrying  Capacity */ /* K_R */ 
      
      Table->Index[n++]   = 11; /* Resource Local Growth Rate */ /* Beta_R */ /* Q */

      Table->Index[n++]   = 12; /* External Immigration Rate */ /* F */ /* Lambda_C_0 */
      
      Table->Index[n++]   = 13; /* Consumer Death Rate */ /* F */ /* Delta_C_0 */    
      Table->Index[n++]   = 15; /* WF Death Rate */       /* WF *//* Delta_C_1 */

      Table->Index[n++]   = 16; /* Consumer Attack Rate */ /* Alpha */
      Table->Index[n++]   = 17; /* Nu = 1/Tau, Tau, Handling Time */ /* Nu */
      
      Table->Index[n++]   = 20; /* Consumer Movement Rate  */ /* mu_C */ /* mu_F */

      Table->Index[n++]   = 29; /* Establishment Rate  (nest starting) *//* W ---> Q */

      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    case 19: /* AZTECA_4D_1 * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every/some of the species: */
      Table->No_of_EVENTS = 3;  /* (Only Diffusion + External Immigration + Death) 
				                            W, Q, and F can undergo these three same processes 
				                        */
      Table->TOTAL_No_of_EVENTS = 3 * Table->No_of_EVENTS + 5;
      Table->LOCAL_STATE_VARIABLES = 4; /* 1 W + 1 Q + 1 F + 1 WF           */
                                        /* WF \equiv RA (handling consumer or
                                           fly larva developing )     
                                        */

      assert(Table->No_of_RESOURCES == 1);
      
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
	    
      /* Conventions */
      Table->K   = n-1;     /* Label last class */
      Table->W = 0; Table->Q = 1;  Table->F = 2; Table->WF = 3; 

      /* List of the 16 (potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;

      Table->Index[n++] = 5; /* No_of_RESOURCES */

      Table->Index[n++] = 0;  /* Worker Movement/Diffusion Rate */ /* Mu */ 
      Table->Index[n++] = 14; /* Queen Movement/Diffusion Rate  */ /* Lambda_C_1 */ 
      Table->Index[n++] = 20; /* Fly Movement/Diffusion RAte    */ /* Mu_C */
      
      Table->Index[n++] = 6;  /* External Immigration Rate */ /* W */ /* Lambda_R_0 */
      Table->Index[n++] = 8;  /* External Immigration Rate */ /* Q */ /* Lambda_R_1 */
      Table->Index[n++] = 12; /* External Immigration Rate */ /* F */ /* Lambda_C_0 */

      Table->Index[n++] = 10; /* Queen Carrying  Capacity */ /* K_R */

      Table->Index[n++] = 7; /* Death Rate (0) */ /* W */   /* Delta_R_0 */
      Table->Index[n++] = 9; /* Death Rate (1) */ /* Q */   /* Delta_R_1 */
      Table->Index[n++] = 13; /* Consumer Death Rate */ /* F */ /* Delta_C_0 */    
      Table->Index[n++] = 15; /* WF Death Rate */       /* WF *//* Delta_C_1 */
       
      Table->Index[n++] = 11; /* Worker Local Production Rate */ /* Beta_R *//* W ---> W + 1 */
      Table->Index[n++] = 29; /* Worker Local Production Rate */ /* Eta_R  *//* Q ---> Q + 1 */

      Table->Index[n++] = 16; /* Consumer Attack Rate */ /* Alpha */
      Table->Index[n++] = 17; /* Nu = 1/Tau, Tau, Handling Time */ /* Nu */
            
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

  
  case 20: /* ECOEVO_PLANTS * * * * * * * * * */

      /* No_of_EVENTS, i.e, common events that can occur to every/some of the species: */
      Table->No_of_EVENTS = 8;  /* All species (strains or phenotypes) can undergo the 
                                   same 8 processes  
				                        */
      Table->TOTAL_No_of_EVENTS = Table->No_of_EVENTS * Table->No_of_RESOURCES;
      Table->LOCAL_STATE_VARIABLES = 2 * Table->No_of_RESOURCES; /* 1 RP + 1 R */
                                        /* RP, resource propagules
                                           R,  individual resource items      
                                        */
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
	    
      /* Conventions */
      Table->K   = n-1; /* Label last class */
      Table->RP = 0; Table->R = 1; 

      /* List of the 16 (potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;

      Table->Index[n++] = 5;                                           /* No_of_RESOURCES */

      Table->Index[n++] = 0;  /* Propagule Movement/Diffusion Rate */  /* Mu */ 
      
      Table->Index[n++] = 6;  /* External Immigration Rate */ /* RP */ /* Lambda_R_0 */
      
      Table->Index[n++] = 10; /* Queen Carrying  Capacity */           /* K_R */

      Table->Index[n++] = 7; /* Death Rate (0) */ /* RP */             /* Delta_R_0 */
      Table->Index[n++] = 9; /* Death Rate (1) */ /* R  */             /* Delta_R_1 */
       
      Table->Index[n++] = 13; /* Propagule Local Death Rate */         /* Delta_C_0 *//* RP ---> RP - 1 */
      Table->Index[n++] = 15; /* Adult Local Death Rate */             /* Delta_C_1 *//* R  ---> R  - 1 */

      Table->Index[n++] = 11; /* Propagule Local Production Rate */    /* Beta_R *//* RP ---> RP + 1 */
      Table->Index[n++] = 29; /* Propagule Local Establisment Rate */  /* Eta_R  *//* R  ---> R + 1 */
      
      Table->Index[n++] = 27; /* Mutation Probability */                 /* p_1 */
      Table->Index[n++] = 28; /* Tradeoff Constant, for instance, R_0 */ /* p_2 */
            
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

  case 21: /* ECO_PLASMIDS * * * * * * * * * */

  #if defined OPLUS
      /* No_of_EVENTS, i.e, common events that can occur to every species: */
      Table->No_of_EVENTS = 7;  /* All species (strains or phenotypes) can undergo the 
                                   same 7 processes  */
      Table->TOTAL_No_of_EVENTS = (Table->No_of_EVENTS-1)*Table->No_of_RESOURCES + Table->No_of_CONJUGATION_EVENTS; 
  #else 
      /* No_of_EVENTS, i.e, common events that can occur to every species: */
      Table->No_of_EVENTS = 7;  /* All species (strains or phenotypes) can undergo the 
                                   same 7 processes  */
      Table->TOTAL_No_of_EVENTS    = Table->No_of_EVENTS * Table->No_of_RESOURCES;
  #endif    
      
      Table->LOCAL_STATE_VARIABLES = Table->No_of_RESOURCES; /* No of different subpopulations */
                                        
      n = 0;
      for(i=0; i<Table->No_of_CELLS; i++)
	      for(j=0; j<Table->LOCAL_STATE_VARIABLES; j++)
	        n++;
	    
      /* Conventions */
      Table->K   = n-1; /* Label last class */
      Table->R   = 0;  

      /* List of the 15 (potentially searcheable) model parameters:  */
      /* MODEL_PARAMETER_SPACE_MAXIMUM (see MODEL_DEFINE_MAX....h) */
      n = 0;
      Table->Index[n++] = 1;  /* No_of_INDIVIDUALS: No_of_PLASMIDS */                                          /* -HN: No_of_PLASMIDS */
      /* No_of_PLASMIDS will be initilized with this argument (No_of_INDIVIDUALS) */

      Table->Index[n++] = 5;  /* No_of_RESOURCES: No_of_SPECIES or bacteria No_of_STRAINS */                   /* -HS:  No_of_STRAINS */
      /* No_of_STRAINS will initialized with this argument (No_of_RESOURCES), but then 
         No_of_RESOURCES will be recalculated according to the contraints 
         given by the interaction structure of the interaction matrices 
         (see: int Determining_actual_No_of_RESOURCES(...)::Patameter_Table.c)
      */
      Table->Index[n++] = 10; /* Carrying  Capacity */ /* K_R */                                               /* -HK:  K_R */

      Table->Index[n++] = 7;  /* Death Rate (0) */ /* Delta_R_0 */                                             /* -H1:  Delta_R_0 */ 
                             
                              /* Delta = Delta_R_0 + Delta_R_1 * (1.0-Resistance ) */

      Table->Index[n++] = 9;  /* Death Rate (1) */ /* Delta_R_1 */                                             /* -H3:  Delta_R_1 */

      Table->Index[n++] = 13; /* Competition induced Death Rate (competition matrix) */                        /* -H6:  Delta_C_0 */      
      
      Table->Index[n++] = 11; /* Cell division rate *//* R ---> R + 1 */                                       /* -H4:  Beta_R */

      Table->Index[n++] = 16; /* Reproduction cost associated to the presence of a plasmid */                  /* -H9:  Alpha_C_0 */ 

      Table->Index[n++] = 17; /* Resistance to stress-induced death associated to the presence of a plasmid */ /* -H10: Nu_C_0    */

      Table->Index[n++] = 0;  /* Bacteria Movement/Diffusion Rate */                                           /* -Hu:  Mu */ 
      
      Table->Index[n++] = 6;  /* External Immigration Rate */                                                  /* -H0:  Lambda_R_0 */
      
      Table->Index[n++] = 8;  /* Conjugation/encounter rate */                                                 /* -H2:  Lambda_R_1 */    
      
      Table->Index[n++] = 18; /* Plasmid Transmission Probability, x */                                        /* -H11:  Xhi_C_0 */
      
      Table->Index[n++] = 27; /* Segregation error at reproduction */                                          /* -Hp1:  p_1 */
       
      Table->Index[n++] = 28; /* Sparsity parameter (a matrix connectance ) */                                 /* -Hp2:  p_2 */ 
            
      Table->TOTAL_No_of_MODEL_PARAMETERS = n;
      break;

    default:
      printf(" This TYPE_of_MODEL (%d) code is not defined.\n", TYPE_of_MODEL);
      printf(" Models (0 to 10): Check input argument list!!!\n");
      exit(0);
   }
   
   Table->TOTAL_GRAND_No_of_EVENTS = Table->TOTAL_No_of_EVENTS * Table->No_of_CELLS;
  
  /* Conventionally, the last label in the argument list of 

     Model_Variables_Code_into_Parameter_Table (...),

     (*K), should be the label of the last model state variable.
     Then ( * K) + 1 becomes de total number of dynamic variables.
  */
  printf("Total No of MODEL PARAMETERS : %d\n", Table->TOTAL_No_of_MODEL_PARAMETERS);
  Print_Press_Key(1,0,"."); 
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

    case 15: /* DIFFUSION_STOLLENBERG_4D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	      for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	        n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 16: /* DIFFUSION_HII_nD * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	      for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	        n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 17: /* DIFFUSION_AZTECA_4D * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	      for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	        n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 18: /* DIFFUSION_AZTECA_4D_0 * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	      for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	        n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 19: /* DIFFUSION_AZTECA_4D_0 * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	      for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	        n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;

    case 20: /* DIFFUSION_ECOEVO_PLANTS * * * * * * * * * * * * * * * * * * * * * * */
      
      n = 0;
      for(i=0; i<P->No_of_CELLS; i++)
	      for(j=0; j<P->LOCAL_STATE_VARIABLES; j++)
	        n++;
       
      /* Conventions */
      P->K   = n-1;     /* Label last class            */
      
      break;
    
    case 21: /* DIFFUSION_ECO_PLASMIDS * * * * * * * * * * * * * * * * * * * * * * */
      
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
