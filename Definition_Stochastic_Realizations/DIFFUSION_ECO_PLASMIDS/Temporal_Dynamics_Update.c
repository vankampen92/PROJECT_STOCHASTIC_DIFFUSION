/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2021 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#define TOLERANCE 1.0E-9


void Updating_Event_Delta_Matrix_Partial_Events(Community * Pa, int Type_of_Event, int Sp, 
                                                Parameter_Table * Table);

/* This functions allocate, initialize and free a number of local communities,
   which make up our total patch system or metapopulation */
extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Temporal_Dynamics_Update( Community ** My_Community,
			                         Parameter_Table * Table,
			                         Stochastic_Rate * Rate,
			                         int Type_of_Event, int * Patch )
{
  /* This function calculates the stochastic rates after the execution of a stochastic event
     in terms of the old ones, with no recalculation. This is a way to optimize the algorithm.
     It is always worth trying this optimization, particularly, for very sparsely coupled 
     spatial systems. 
  */
  /* Input arguments:
     . My_Community is a pointer to the whole patch system
     . Table        is a pointer to a Parameter_Table type of structure from which other 
                    structures hang. 
     . Rate         is a pointer to the structure 'Stocastic Rate' 
     . Type_of_Event is an index to the local event occurring in a Patch 
                    (between 0 and Table->TOTAL_No_of_EVENTS-1)
     . Patch        is an array containing the two patches involved in a movement event. 
                    Patch[0] is the patch sending the individual
		                Patch[1] is the patch receiving the individual.
                    Patch[2] is the strain ID of the RECIPIENT species
                    Patch[3] is the strain ID of the TRANSCONJUGANT species
                    Patch[4] is 0 if no update is required and 1 otherwise. 

     Output arguments: 
     . Rate         Stochastic Rate is updated from previous value (without recalculating)
  */
  int Type_Event;      /* Index of the event between 0 and Table->No_of_EVENTS-1   */
                       /* Type of Event that occurs at a given species on patch x  */
  bool bool_Next_Time;
  int i, n, m, k; 
  int x, y, Sp, Sp_0, Sp_Recip, Sp_Trans; 
   
  double Delta_Rate, time_old, lambda_old, lambda_new, time_current, Next_Time;  
  double OutMigration; 

  assert(Patch[4] == 1);    /* Update of temporal update required !!! */

  n = Table->No_of_EVENTS;  /* Number of Event a given Species can undergo (0 to 6) */

  /* Remember: 
              TOTAL_No_of_EVENTS = (No_of_EVENTS-1)*No_of_RESOURCES + No_of_CONJUGATION_EVENTS; 
  */
  x = Patch[0]; 
  y = Patch[1]; 
  
  Community * Pa = My_Community[x];

  if ( Type_of_Event >= Table->TOTAL_No_of_EVENTS ) {
    printf(" Error in Temporal Dynamics Update!!!\n");
    printf(" Type of Event occurring is too large.\n");
    printf(" It can only be labeled from 0 or %d\n", Table->TOTAL_No_of_EVENTS-1);
    printf(" but events is seemingly labeled as %d\n", Type_of_Event); 
    printf(" The program will exit\n"); 
    exit(0); 
  }
 
  if (Patch[4] == 1) { /* If Patch[4] is 0, no update of the array controling total rates 
                          is required */

    if ( x == y ) {  /* LOCAL PROCESS has happened on a single Patch x */
      if ( Type_of_Event >= Table->No_of_RESOURCES * (Table->No_of_EVENTS - 1) ) {
        /* (Event is 6: Conjugation Event) */
        Sp_Recip = Patch[2];
        Sp_Trans = Patch[3];
        assert( Sp_Trans != Sp_Recip);
        Delta_Rate = 0.0; bool_Next_Time = false;
        Temporal_Dynamics_Update_Conjugation( Pa, x, Table, Rate, Type_of_Event, Patch,
                                              &bool_Next_Time, &Delta_Rate );
      }
      else {
        /* (Event can only be from 1 to 5) */ 
        Sp         = Type_of_Event/(Table->No_of_EVENTS-1);
        Type_Event = Type_of_Event%(Table->No_of_EVENTS-1);
  
        assert( Type_Event > 0 && Type_Event < 6 );  /* Because this Type_Event = 0 is an out
							                                          migration events. Only
							                                          possible when x patch 
							                                          is different from y patch, 
                                                        and Type_Event = 6 is a conjugation event */                                
      if(Type_Event == 5) {
      /* Sp_0 is the Strain ID of the plasmid free profile strain of the same bacterial type */
        Sp_1 = Patch[2];
        k    = Table->StrainType_and_Profile[Sp_1][0];
        Sp_0 = Table->n_0[k]; 
        assert(Sp_1 == Sp_0);
        Sp  = Sp_1;          /* The Strain ID of a plasmid free bacterial type of type k */ 
        Type_Event = 4; 
      }

      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Focal_Species ( Pa, x, Table, Rate, Type_of_Event,  
                                               Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing conjugation rates of the focal species with the species in its recipient list, 
         and the conjugation rates of each of species in the donor list of the focal species 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Conjugation_Rates_Focal_Species ( Pa, x, Table, Rate, Type_of_Event,  
                                                                 Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing rates across species/types as a result 
         of a change in the amount of empty space (m_0). 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Empty_Space_Induced ( Pa, x, Table, Rate, Type_of_Event,
                                                     Type_Event, bool_Next_Time, &Delta_Rate );
      /* Changing rates of the competition-induced death of the species competing with
         the focal one (in ist competition list).  Updating Per Capita Competition rates. 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Competition_Rates ( Pa, x, Table, Rate, Type_of_Event,
                                                   Type_Event, bool_Next_Time, &Delta_Rate );
      } 
    }
    else {  /* MOVEMENT EVENT involving two Patches: 
             x: patch exporting an individual propagule 
             y: patch receiving an individual propagule, 
             out from patch 'x' into patch 'y'   
          */       
		      /* Out migration event sending one individual:
             RP: propagules
          */
      Sp         = Type_of_Event/(Table->No_of_EVENTS-1);
      Type_Event = Type_of_Event%(Table->No_of_EVENTS-1);         

      assert( Type_Event == 0 ) ; 
      assert( Table->No_of_CELLS > 1 );

      /* Changes on Patch x: 
      */    
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Focal_Species ( Pa, x, Table, Rate, Type_of_Event,  
                                             Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing conjugation rates of the focal species with those in it recipient list, 
         and the conjugation rates of each of species in the donor list of the focal one 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Conjugation_Rates_Focal_Species ( Pa, x, Table, Rate, Type_of_Event,  
                                                               Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing rates of all other species/types as a result of a change in the amount 
         of empty space (m_0) induced by the evant that has occurred. 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Empty_Space_Induced ( Pa, x, Table, Rate, Type_of_Event,
                                                   Type_Event, bool_Next_Time, &Delta_Rate );
      /* Changing rates of the competition-induced death of the species competing with
         the focal one (in ist competition list).  Updating Per Capita Competition rates. 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Competition_Rates ( Pa, x, Table, Rate, Type_of_Event,
                                                 Type_Event, bool_Next_Time, &Delta_Rate );    
    /* Changes on Patch y: 
    (in the rates of all the events due to the adquisition of an individual of species Sp in patch y 
    */
      Pa         = My_Community[y];
      Type_Event = 1; /* Immigration into patch y */

      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Focal_Species ( Pa, y, Table, Rate, Type_of_Event,  
                                               Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing conjugation rates of the focal species with those in it recipient list, 
         and the conjugation rates of each of species in the donor list of the focal one 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Conjugation_Rates_Focal_Species ( Pa, y, Table, Rate, Type_of_Event,  
                                                               Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing rates of all other species/types as a result 
         of a change in the amount of empty space (m_0). 
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Empty_Space_Induced ( Pa, y, Table, Rate, Type_of_Event,
                                                   Type_Event, &bool_Next_Time, &Delta_Rate );
      /* Changing rates of the competition-induced death of the species competing with
         the focal one (in ist competition list). Updating Per Capita Competition rates.  
      */
      Delta_Rate = 0.0; bool_Next_Time = false;
      Temporal_Dynamics_Update_Competition_Rates ( Pa, y, Table, Rate, Type_of_Event,
                                                 Type_Event, &bool_Next_Time, &Delta_Rate );
    }
  }
  else {
    assert ( Type_of_Event >= Table->No_of_RESOURCES * (Table->No_of_EVENTS - 1));

    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        /* The time-to-happen of the conjugation that failed is recalculted 
           and reordered ('bubbled') throughout the priority quee binary tree 
        */
        int Index_Node_x = Table->TOTAL_No_of_EVENTS * x + Type_of_Event;
        treenode * Node  = Table->Tree_Node_Index[Index_Node_x];
        Node->value = time_current - 1.0/Pa->rToI[Type_of_Event] * log(RANDOM); 
        bubbling(Node, Table->Tree_Node_Index);
    #endif

    printf("Unsuccessful Conjugation Event...\n"); 
  } 

  if(Rate->Total_Rate <= 0.0){
    printf("\n");
    printf(" R is the total temporal rate of system configuration change\n");
    printf(" R = %g\n", Rate->Total_Rate );
    printf(" If R is zero, no further change is possible\n");
    printf(" but R shouldn't be too negative!!!\n");
    printf(" If it is, check if it is lower than certain very small tolerance value\n");
    printf("\n");
    // if( Rate->Total_Rate < -TOLERANCE) exit(0);
  }
}

void Updating_Event_Delta_Matrix_Partial_Events(Community * Pa, int Type_of_Event, int Sp, 
                                                Parameter_Table * Table)
{
  /* 
     This is the subset of the Delta Matrix entries that depend
     on current system configuration for the first 5 events (from 0 to 4). 
     Therefore, they need to be changed in agreement to the event that 
     has just occurred. 
  */
  int i, n_Event; 

  double ** Delta_Matrix = Pa->Event_Delta_Tensor[Sp];

  int * n = Pa->n;

  double K_R = Table->K_R;
  int m_0    = Pa->m_0;
  double B0  = Table->Beta_AP[Sp] * (1.0 - Table->p_1);
  double B1  = Table->Beta_AP[Sp] * Table->p_1;

  n_Event = Sp * (Table->No_of_EVENTS-1) + 3;  /* Type_of_Event = 3: Competition-Induced Death */ 

  assert(n_Event < Table->No_of_RESOURCES * (Table->No_of_EVENTS-1) );

  double C_Alpha = Pa->rate[n_Event];  /* Per capita Rate */

  assert( Table->R == 0 ); 
  
  switch( Type_of_Event )
    {
    case 0:  /* 0: Out-Migration (R --> R-1) and some other patch gains one */
       
      Delta_Matrix[0][3] = -C_Alpha; 
      Delta_Matrix[0][4] =  B0/K_R * ((double)n[Sp] - (double)m_0 + 1.0); 
      Delta_Matrix[0][5] =  B1/K_R * ((double)n[Sp] - (double)m_0 + 1.0);
    break;
      
    case 1:  /* 1: External Immigration event */
      Delta_Matrix[1][3] = +C_Alpha; 
      Delta_Matrix[1][4] =  B0/K_R * ((double)m_0 - (double)n[Sp] + 1.0); 
      Delta_Matrix[1][5] =  B1/K_R * ((double)m_0 - (double)n[Sp] + 1.0);
      
    break;
    
    case 2: /* 2: Death  */
      Delta_Matrix[2][3] = -C_Alpha; 
      Delta_Matrix[2][4] =  B0/K_R * ((double)n[Sp] - (double)m_0 + 1.0); 
      Delta_Matrix[2][5] =  B1/K_R * ((double)n[Sp] - (double)m_0 + 1.0);
      
    break;
    
    case 3: /* 3: Competition-Induced Death */
      Delta_Matrix[3][3] = -C_Alpha; 
      Delta_Matrix[3][4] =  B0/K_R * ((double)n[Sp] - (double)m_0 + 1.0); 
      Delta_Matrix[3][5] =  B1/K_R * ((double)n[Sp] - (double)m_0 + 1.0);
      
    break;
      
    case 4: /* 4: Cell Division (withouth error)*/
      Delta_Matrix[4][3] = +C_Alpha; 
      Delta_Matrix[4][4] =  B0/K_R * ((double)m_0 - (double)n[Sp] + 1.0); 
      Delta_Matrix[4][5] =  B1/K_R * ((double)m_0 - (double)n[Sp] + 1.0);
    
    break;
        
    default:
      /* Something is very very wrong!!! */
      printf(" Type_of_Event = %d\t This value is not possible!!!\n", Type_of_Event);
      printf(" Only 0 to 4 are possible. Type of Event is ill-defined\n");
      printf(" The program will exit\n");
      Print_Press_Key(1,0,"."); 
      exit(0);
    }
}
