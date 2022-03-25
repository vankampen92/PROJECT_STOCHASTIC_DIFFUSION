/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2021 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#define TOLERANCE -1.0E-10

void Updating_Event_Delta_Matrix(Community * Pa, int Type_of_Event, Parameter_Table * Table);

void Temporal_Dynamics_Update( Community ** My_Community,
			       Parameter_Table * Table,
			       Stochastic_Rate * Rate,
			       int Type_of_Event, int * Patch)
{
  /* This function calculates the stochastic rates after the execution of a stochastic event
     in terms of the old ones, with no recalculation. This is a way to optimize the algorithm.
     It is always worth trying this optimization, particularly, for very sparsely coupled 
     systems. 
  */
  /* Input arguments:
     
     . My_Community is a pointer to the whole patch system
     . Table        is a pointer to a Parameter_Table type of structure from which other 
                    structures hang. 
     . Rate         is a pointer to the structure Stocastic Rate 
     . Type_of_Event is a label to the event occurring in a Patch 
     . Patch        is an array containing the two patches involved in a movement event. 
                    Patch[0] is the patch sending the individual
		    Patch[1] is the patch receiving the individual. 

     Output arguments: 
     . Rate         Stochastic Rate is updated from previous value (without recalculating)
  */
  int i, n, m, k; 
  int x, y; 
  double Delta_Rate;
  
  Community * Pa;
  double OutMigration; 

  n = Table->TOTAL_No_of_EVENTS;
  
  x = Patch[0]; y = Patch[1];
  
  if ( x == y ) { /* LOCAL PROCESS on a single Patch x */

    if (Type_of_Event < n ) {
      
      assert( Type_of_Event != 0 && Type_of_Event != 3 ) ; /* Because these are out
							      migration events. Only
							      possible when x patch 
							      is different from y patch
							   */
      
      Pa    = My_Community[x];
      Updating_Event_Delta_Matrix(Pa, Type_of_Event, Table); 
      
      m = Pa->Event_Adjacence_List[Type_of_Event][n]; /* How many events are connected 
							 to event 'Type_of_Event'???
							 i.e., the lenght of the adjacence  
							 list of 'Type_of_Event'
						      */
      Delta_Rate = 0.0; 
      for(i=0; i<m; i++) {
	k = Pa->Event_Adjacence_List[Type_of_Event][i]; /* Which events are connected 
							   to event 'Type_of_Event'???
							*/
	Delta_Rate  += Pa->Event_Delta_Matrix[Type_of_Event][k];
	Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event][k];
      }

      Pa->ratePatch    += Delta_Rate; 				
      Rate->Total_Rate += Delta_Rate; 				
      Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch );
    }
    else {
      printf(" Error in Temporal Dynamics Update!!!\n");
      printf(" Type of Event occurring is too large.\n");
      printf(" It can only be labeled from 0 or %d\n", m-1);
      printf(" but events is seemingly labeled as %d\n", Type_of_Event); 
      printf(" The program will exit\n"); 
      exit(0); 
    }   
  }
  else {         /* MOVEMENT EVENT involving two Patches */
		 /* Out migration sending one individual (R or C, 
		    RA individuals do not move) out 
		    from patch 'x' into patch 'y'      */
    
    assert( Type_of_Event == 0 || Type_of_Event == 3 ) ; 
    
    /* Changes in rates due to the loss of an individual in patch x */
    Pa    = My_Community[x];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event, Table);
    
    m = Pa->Event_Adjacence_List[Type_of_Event][n]; /* How many events are connected 
						       to event 'Type_of_Event'???
						       Lenght of the adjacence list 
						       of 'Type_of_Event'
						    */
    Delta_Rate = 0.0;
    for(i=0; i<m; i++) {
      k = Pa->Event_Adjacence_List[Type_of_Event][i]; /* Which events are connected 
							 to event 'Type_of_Event'???
						      */
      Delta_Rate  += Pa->Event_Delta_Matrix[Type_of_Event][k];
      Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event][k];
    }
    Pa->ratePatch    += Delta_Rate; 				
    Rate->Total_Rate += Delta_Rate; 				
    Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch );
      
    
    /* Changes in rates due to the adquisition of an individual in patch y */
    Pa    = My_Community[y];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event+1, Table);
    
    m = Pa->Event_Adjacence_List[Type_of_Event+1][n]; /* How many events are connected 
							 to event 'Type_of_Event+1'???
							 Lenght of the adjacence list 
							 of 'Type_of_Event+1', 
							 because 1 or 4 are events  
							 involving the increase of 
							 the resource or the consumer,
							 respectively
						    */
    Delta_Rate = 0.0;
    for(i=0; i<m; i++) {
      k = Pa->Event_Adjacence_List[Type_of_Event+1][i]; /* Which events are connected 
							   to event 'Type_of_Event'???
							*/
      Delta_Rate  += Pa->Event_Delta_Matrix[Type_of_Event+1][k];
      Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event+1][k];
    }
    Pa->ratePatch    += Delta_Rate; 				
    Rate->Total_Rate += Delta_Rate; 				
    Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch );
  }
  
  if(Rate->Total_Rate <= 0.0){
    printf("\n");
    printf(" R is the total temporal rate of system configuration change\n");
    printf(" R = %g\n", Rate->Total_Rate );
    printf(" If R is zero, no further change is possible\n");
    printf(" but R shouldn't be too negative!!!\n");
    printf(" If it is, check if it is lower than certain very small tolerance value\n");
    printf("\n");
    // if( Rate->Total_Rate < TOLERANCE) exit(0);
  }
}

void Updating_Event_Delta_Matrix(Community * Pa, int Type_of_Event, Parameter_Table * Table)
{
  /* 
     This is the subset of the Delta Matrix entries that depend
     on system configuration. Therefore, they need to be 
     changed in agreement to the event that has just occurred 
  */
  
  
  double ** Delta_Matrix = Pa->Event_Delta_Matrix;

  int * n = Pa->n;

  double K_R = Table->K_R;

  /* Definition of the local state vector numerical order */
  #include <Model_Variables_Code.Include.c>
  
  switch( Type_of_Event )
    {
    case 0:  /* Resource Out-Migration (R --> R-1) and some other patch gains one */ 
      Delta_Matrix[0][6] = -Table->Alpha_C_0/K_R * (double)n[A];
      Delta_Matrix[0][7] = Table->Beta_R/K_R * (2.0*(double)n[R]-K_R+1.0); 
    break;
      
    case 1:  /* Resource External Immigration event */
      Delta_Matrix[1][6] = Table->Alpha_C_0/K_R * (double)n[A];
      Delta_Matrix[1][7] = Table->Beta_R/K_R * (K_R-2.0*(double)n[R]+1.0); 
    break;
    
    case 2:  /* Resoure Death  */
      Delta_Matrix[2][6] = -Table->Alpha_C_0/K_R * (double)n[A];;
      Delta_Matrix[2][7] = Table->Beta_R/K_R * (2.0*(double)n[R]-K_R+1.0); 
    break;
    
    case 3:  /* Consumer Out-Migration (A --> A-1) and some other patch gains one */ 
      Delta_Matrix[3][6] = -Table->Alpha_C_0/K_R * (double)n[R];
      
    break;
      
    case 4:  /* Consumer External Immigration event  */  
      Delta_Matrix[4][6] = Table->Alpha_C_0/K_R * (double)n[R];
      
    break;
    
    case 5:  /* Consumer Death  */
      Delta_Matrix[5][6] = -Table->Alpha_C_0/K_R * (double)n[R];
      
    break;
    
    case 6:  /* Consumer Consumption of resource and dimmer formation */
      Delta_Matrix[6][6] = -Table->Alpha_C_0/K_R * (1.0+(double)n[A]+(double)n[R]);
      Delta_Matrix[6][7] = Table->Beta_R/K_R * (2.0*(double)n[R]-K_R+1.0); 
    break;

    case 7:  /* Local Growthh of Resources */
      Delta_Matrix[7][6] = Table->Alpha_C_0/K_R * (double)n[A];
      Delta_Matrix[7][7] = Table->Beta_R/K_R * (K_R-2.0*(double)n[R]+1.0); 
    break;

    case 8:  /* Local Growth of Consumers  */
      Delta_Matrix[8][6] = Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    case 9:  /* RA Local Death: RA ---> RA - 1 */
             /* Delta Matrix change is independent from current system configuration */
    break;

    case 10: /* RA relax back into A: RA ---> A */
      Delta_Matrix[10][6] = Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    default:
      /* Something is very very wrong!!! */
      printf(" Type_of_Event = %d\t This value is not possible!!!\n", Type_of_Event);
      printf(" Type of Event ill-defined\n");
      printf(" The program will exit\n");
      Press_Key(); 
      exit(0);
    }
}
