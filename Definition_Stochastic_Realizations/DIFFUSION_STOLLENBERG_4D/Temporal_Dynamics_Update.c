/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2021 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#define TOLERANCE -1.0E-10

void Updating_Event_Delta_Matrix(Community * Pa, int Type_of_Event, Parameter_Table * Table);

void Temporal_Dynamics_Update( Community ** My_Community,
			                         Parameter_Table * Table,
			                         Stochastic_Rate * Rate,
			                         int Type_of_Event, int * Patch )
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
  treenode * Leaf; 
  
  Community * Pa;
  double OutMigration; 

  n = Table->TOTAL_No_of_EVENTS;
  
  x = Patch[0]; y = Patch[1];
  
  if ( x == y ) {  /* LOCAL PROCESS on a single Patch x */
    if (Type_of_Event < n ) {
      
      assert( Type_of_Event != 0 && Type_of_Event != 6 ) ; /* Because these are out
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
      
      #if defined BINARY_TREE_OPTIMIZATION
        /* Updating of the Binary Tree to Sample Discrete Distribution */
        Leaf = Table->Leaves[x];
        sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Rate);
      #endif
      
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
  else {  /* MOVEMENT EVENT involving two Patches */       
		      /* Out migration sending one individual (RP or C, 
		         RA individuals do not move) out 
		         from patch 'x' into patch 'y'      
          */
    assert( Type_of_Event == 0 || Type_of_Event == 6 ) ; 
    
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
    
    #if defined BINARY_TREE_OPTIMIZATION 
      /* Updating of the Binary Tree to Sample Discrete Distribution */  
      Leaf = Table->Leaves[x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Rate);  
    #endif

    Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch ); 

    /* Changes in rates due to the adquisition of an individual in patch y */
    Pa    = My_Community[y];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event+1, Table);
    
    m = Pa->Event_Adjacence_List[Type_of_Event+1][n]; /* How many events are connected 
							                                           to event 'Type_of_Event+1'???
							                                           Lenght of the adjacence list 
							                                           of 'Type_of_Event+1', 
							                                           because 1 or 7 are events  
							                                           involving the increase of 
							                                           the propagule or the consumer
							                                           population, respectively
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

    #if defined BINARY_TREE_OPTIMIZATION
      /* Updating of the Binary Tree to Sample Discrete Distribution */
      Leaf = Table->Leaves[y];  
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Rate);
    #endif

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
    case 0:  /* 0: Propagule Out-Migration (P --> P-1) and some other patch gains one */ 
      Delta_Matrix[0][4] = -Table->Eta_R/K_R * (K_R - (double)n[R]);
      
    break;
      
    case 1:  /* 1: Propagule External Immigration event */
      Delta_Matrix[1][4] = +Table->Eta_R/K_R * (K_R - (double)n[R]);
      
    break;
    
    case 2: /* 2: Propagule Death  */
      Delta_Matrix[2][4] = -Table->Eta_R/K_R * (K_R - (double)n[R]);
      
    break;
    
    case 3: /* 3: Propagule Production */
      Delta_Matrix[3][4] = +Table->Eta_R/K_R * (K_R - (double)n[R]);
      
    break;
      
    case 4: /* 4: Propagule Establishment  */  
      Delta_Matrix[4][4] = -Table->Eta_R/K_R * (1.0 + (double)n[RP] + K_R - (double)n[R]);
      Delta_Matrix[4][9] = +Table->Alpha_C_0/K_R * (double)n[A];
      
    break;
    
    case 5:  /* 5: Resource Death  */
      Delta_Matrix[5][4] = +Table->Eta_R/K_R * (double)n[RP];
      Delta_Matrix[5][9] = -Table->Alpha_C_0/K_R * (double)n[A];
       
    break;
    
    case 6:  /* 6: Consumer Out-Migration (A ---> A-1) and some other patch gains one */ 
      Delta_Matrix[6][9] = -Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    case 7:  /* 7: Consumer External Immigration (A ---> A + 1) */  
      Delta_Matrix[7][9] = +Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    case 8:  /* 8: Searching Consumer Death      (A ---> A - 1)*/   
      Delta_Matrix[8][9] = -Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    case 9:  /* 9: Attack: Consumer Consumption of a Resource (A + R ---> RA) */
      Delta_Matrix[9][4] = +Table->Eta_R/K_R * (double)n[RP];
      Delta_Matrix[9][9] = -Table->Alpha_C_0/K_R * (1.0 + (double)n[R] + (double)n[A]);
      
    break;

    case 10: /* 10: Handling Consumers relax back into Free Consumers ( RA ---> A) */ 
      Delta_Matrix[10][9] = Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    case 11: /* 11: Consumer Growth:  (RA ---> RA + A)  */
      Delta_Matrix[11][9] = Table->Alpha_C_0/K_R * (double)n[R];
      
    break;

    case 12: /* 12: (Handling) Consumer Death (RA ---> RA-1) */
      /* Delta values are independent from system configuration */
      
    break;

    default:
      /* Something is very very wrong!!! */
      printf(" Type_of_Event = %d\t This value is not possible!!!\n", Type_of_Event);
      printf(" Type of Event ill-defined\n");
      printf(" The program will exit\n");
      Print_Press_Key(1,0,"."); 
      exit(0);
    }
}
