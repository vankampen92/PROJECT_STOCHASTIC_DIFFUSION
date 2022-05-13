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
     It is always worth trying this optimization for very sparsely coupled systems. 
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
      
      assert( Type_of_Event != 0 ) ; /* Because 0 labels an out
					migration events. Only
					possible when x patch 
					is different from y patch
				     */
      
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
		 /* Out migration sending one individual (C) 
		    out from patch 'x' to patch 'y'      
		 */
    assert( Type_of_Event == 0 ) ; 
    
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
							 because 1 is the event  
							 involving the increase of 
							 the consumer 
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
  double ** Delta_Matrix = Pa->Event_Delta_Matrix;

  int * n = Pa->n;

  double K_R = Table->K_R;

  /* Definition of the local state vector numerical order */
  #include <Model_Variables_Code.Include.c>
  
  switch( Type_of_Event )
    {
    case 0:  /* Consumer outmigration  */ 
      
      Delta_Matrix[0][4] = -Table->Chi_C_0/K_R * (double)n[RA];
      
    break;
      
    case 1:  /* Consumer immigration */

      Delta_Matrix[1][4] = +Table->Chi_C_0/K_R * (double)n[RA];
      
    break;
    
    case 2:  /* Consumer Attack  */
      
      Delta_Matrix[2][4] = Table->Chi_C_0/K_R * (1.0+(double)n[A]-(double)n[RA]);
      
    break;
    
    case 3:  /* Consumer Handling */ 
      
      Delta_Matrix[3][4] = Table->Chi_C_0/K_R * (1.0-(double)n[A]+(double)n[RA]);
      
    break;

    case 4:  /* Consumer Interference (A + RA --> ARA): Triplet formation  */ 
      
      Delta_Matrix[4][4] = -Table->Chi_C_0/K_R * (1.0+(double)n[A]+(double)n[RA]);
      
    break;

    case 5:  /* Consumer Interference (ARA --> A + RA): Triple degration */ 
      
      Delta_Matrix[5][4] = -Table->Chi_C_0/K_R * (1.0-(double)n[A]-(double)n[RA]);
      
    break;
      
    default:
      /* Something is very very wrong!!! */
      printf(" Type_of_Event = %d\t This value is not possible!!!\n", Type_of_Event);
      printf(" There are only six events possible labelled as (0, 1, 2, 3, 4, 5)\n");
      printf(" Type of Event ill-defined\n");
      printf(" The program will exit\n");
      Press_Key(); 
      exit(0);
    }
}
