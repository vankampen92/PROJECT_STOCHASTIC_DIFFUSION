/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2021 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#define TOLERANCE 1.0E-9

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
     . Patch        is an array containing the two patches involved in a movement event. 
                    Patch[0] is the patch sending the individual
		    Patch[1] is the patch receiving the individual. 
     Output arguments: 
     . Rate         Stochastic Rate is updated from previous value (without recalculating)
  */
  bool bool_Next_Time;
  int i, n, m, k; 
  int x, y; 
  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double Delta_Rate;
  double time_old, lambda_old, lambda_new, time_current, Next_Time;
  treenode * Leaf;
  treenode * Node;  
  
  Community * Pa;
  double OutMigration; 

  n = Table->TOTAL_No_of_EVENTS;
  
  x = Patch[0]; y = Patch[1];

  if (Type_of_Event >= n ) {
    printf(" Error in Temporal Dynamics Update!!!\n");
    printf(" Type of Event occurring is too large.\n");
    printf(" It can only be labeled from 0 or %d\n", Table->TOTAL_No_of_EVENTS-1);
    printf(" but events is seemingly labeled as %d\n", Type_of_Event); 
    printf(" The program will exit\n"); 
    exit(0); 
  }

  #if defined BINARY_TREE_SUPER_OPTIMIZATION
    /* The index of the tip the tree corresponding to a particular 
       event in a given cell 
    */
    Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Type_of_Event;  
  #endif 

  if ( x == y ) {  /* LOCAL PROCESS on a single Patch x */
  
    assert( Type_of_Event != 0 && Type_of_Event != 5 ) ; /* Because (0) and (5) are 
                                                            out migration events. Only
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
    bool_Next_Time = false; 
    time_current = Rate->Stochastic_Time; 
    for(i=0; i<m; i++) {
	    k = Pa->Event_Adjacence_List[Type_of_Event][i]; /* ith event connected 
							                                           to event 'Type_of_Event'???
							                                        */
	    Delta_Rate  += Pa->Event_Delta_Matrix[Type_of_Event][k];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + k;
        Leaf = Table->Leaves[Index_Leaf_x];
        sum_Delta_upto_Root(Table->Treeroot, Leaf, 
                            Pa->Event_Delta_Matrix[Type_of_Event][k]);
      #endif                            

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Index_Node_x = Table->TOTAL_No_of_EVENTS * x + k;

        lambda_old   = Pa->rToI[k]; 
        time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
              
        if(time_old == INFINITY) assert(lambda_old == 0.0);
      #endif 

	    Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event][k];

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Node = Table->Tree_Node_Index[Index_Node_x];
        lambda_new   = Pa->rToI[k];
        if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
          Pa->rToI[k] = 0.0;
          lambda_new  = 0.0; 
        }
        if(k == Type_of_Event) {
          if( Pa->rToI[k] == 0.0 ) 
            Next_Time = INFINITY;
          else
            Next_Time = time_current - 1.0/Pa->rToI[k] * log(RANDOM);
                                
          bool_Next_Time = true; 
        }
        else if (lambda_old == 0.0  && lambda_new == 0.0) 
          Next_Time = INFINITY;
        else if (lambda_old > 0.0   && lambda_new == 0.0)
          Next_Time = INFINITY;
        else if (lambda_old == 0.0  && lambda_new > 0.0)
          Next_Time = time_current - 1.0/Pa->rToI[k] * log(RANDOM);
        else if (lambda_old > 0.0  && lambda_new > 0.0) {
          if( time_old <= time_current ) {
            printf("time_old = %g\t time_current = %g\n", time_old, time_current);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) in Patches (%d, %d)\n", 
                    time_current, k, Type_of_Event, Patch[0], Patch[1]); 
            Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
            printtree(Table->Treeroot);

            assert( time_old > time_current );
          }
          Next_Time = time_current + lambda_old/lambda_new * (time_old-time_current);
        }
        else { 
          if(lambda_new < 0.0 || lambda_old < 0.0) {
            printf("A rate has become too negative: lambda_new = %g\t lambda_old = %g\n", 
                    lambda_new, lambda_old);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) in Patches (%d, %d)\n", 
                     time_current, k, Type_of_Event, Patch[0], Patch[1]);
            Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
            printtree(Table->Treeroot);
            Print_Press_Key(1,1,"Kill the program is lambda_new is too negative");   
          }
        }
        Node->value = Next_Time; 
        bubbling(Node, Table->Tree_Node_Index);

        // Print_Press_Key(1,1,"Printing out Tree before after bubbling\n");
        // printtree(Table->Treeroot);
      #endif  
    }

    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      if( bool_Next_Time == false) {
        Index_Node_x = Table->TOTAL_No_of_EVENTS * x + Type_of_Event;
        Node = Table->Tree_Node_Index[Index_Node_x];
        Next_Time = time_current - 1.0/Pa->rToI[Type_of_Event] * log(RANDOM);
        Node->value = Next_Time;         
        bubbling(Node, Table->Tree_Node_Index);
      }
    #endif

    Pa->ratePatch    += Delta_Rate; 				
    Rate->Total_Rate += Delta_Rate;
      
    #if defined BINARY_TREE_OPTIMIZATION
      /* Updating of the Binary Tree to Sample Discrete Distribution */
      Leaf = Table->Leaves[x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Rate);
    #endif
      
    Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch );   
  }
  else {  /* MOVEMENT EVENT involving two Patches: 
             x: patch exporting an individual 
             y: patch receiving an individual, 
             out from patch 'x' into patch 'y'   
          */       
		      /* Out migration event sending one individual:
             0: Plasmid Free Cells at rate Table->Mu
             1: Plasmid Carrying Cells at rate Table->Mu_C
          */  
		      
    assert( Type_of_Event == 0 || Type_of_Event == 5 ) ; 
    assert( Table->No_of_CELLS > 1 );
    
    /* Changes in rates due to 
       the loss of an individual in patch x */
    Pa    = My_Community[x];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event, Table);
    
    m = Pa->Event_Adjacence_List[Type_of_Event][n]; /* How many events are connected 
						                                           to event 'Type_of_Event'???
						                                           m is the lenght of the adjacence list 
						                                           of 'Type_of_Event' */
    Delta_Rate = 0.0;
    bool_Next_Time = false;
    time_current = Rate->Stochastic_Time;
    for(i=0; i<m; i++) {
      k = Pa->Event_Adjacence_List[Type_of_Event][i]; /* Which events are connected 
							                                           to event 'Type_of_Event'???
						                                          */
      Delta_Rate  += Pa->Event_Delta_Matrix[Type_of_Event][k];
      
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + k;
        Leaf = Table->Leaves[Index_Leaf_x];
        sum_Delta_upto_Root(Table->Treeroot, Leaf, 
                            Pa->Event_Delta_Matrix[Type_of_Event][k]);
      #endif

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Index_Node_x = Table->TOTAL_No_of_EVENTS * x + k;

        lambda_old   = Pa->rToI[k]; 
        time_old     = Table->Tree_Node_Index[Index_Node_x]->value;

        if(time_old == INFINITY) assert(lambda_old == 0.0);
      #endif 

      Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event][k];

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Node = Table->Tree_Node_Index[Index_Node_x];
        lambda_new   = Pa->rToI[k];
        if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
          Pa->rToI[k] = 0.0;
          lambda_new  = 0.0; 
        }

        if(k == Type_of_Event) {
          if( Pa->rToI[k] == 0.0 ) 
            Next_Time = INFINITY;
          else
            Next_Time = time_current - 1.0/Pa->rToI[k] * log(RANDOM);
                                
          bool_Next_Time = true; 
        }
        else if (lambda_old == 0.0  && lambda_new == 0.0) 
          Next_Time = INFINITY;
        else if (lambda_old > 0.0   && lambda_new == 0.0)
          Next_Time = INFINITY;
        else if (lambda_old == 0.0  && lambda_new > 0.0)
          Next_Time = time_current - 1.0/Pa->rToI[k] * log(RANDOM);
        else if (lambda_old > 0.0  && lambda_new > 0.0) {
          if( time_old <= time_current ) {
            printf("time_old = %g\t time_current = %g\n", time_old, time_current);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) in Patches (%d, %d)\n", 
                    time_current, k, Type_of_Event, Patch[0], Patch[1]); 
            Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
            printtree(Table->Treeroot);

            assert( time_old > time_current );
          }
          Next_Time = time_current + lambda_old/lambda_new * (time_old-time_current);
        }
        else { 
          if(lambda_new < 0.0 || lambda_old < 0.0) {
            printf("A rate has become too negative: lambda_new = %g\t lambda_old = %g\n", 
                    lambda_new, lambda_old);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) in Patches (%d, %d)\n", 
                     time_current, k, Type_of_Event, Patch[0], Patch[1]);
            Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
            printtree(Table->Treeroot);
            Print_Press_Key(1,1,"Kill the program is lambda_new is too negative");
          }
        }

        Node->value = Next_Time; 
        bubbling(Node, Table->Tree_Node_Index);
      #endif      
    }
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      if( bool_Next_Time == false) {
        Index_Node_x = Table->TOTAL_No_of_EVENTS * x + Type_of_Event;
        Node = Table->Tree_Node_Index[Index_Node_x];
        Next_Time = time_current - 1.0/Pa->rToI[Type_of_Event] * log(RANDOM);
        Node->value = Next_Time; 
        bubbling(Node, Table->Tree_Node_Index);
      }
    #endif 

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
							                                           m is the lenght of the adjacence list 
							                                           of 'Type_of_Event+1', 
							                                           because 1 or 6 are events  
							                                           involving the increase of 
							                                           the plasmid free cell population or 
                                                         the plasmid carrying cell, respectively
						                                          */
    Delta_Rate = 0.0;
    for(i=0; i<m; i++) {
      k = Pa->Event_Adjacence_List[Type_of_Event+1][i]; /* Which events are connected 
							                                             to event 'Type_of_Event'???
							                                          */
      Delta_Rate  += Pa->Event_Delta_Matrix[Type_of_Event+1][k];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Index_Leaf_y = Table->TOTAL_No_of_EVENTS * y + k;  
        Leaf = Table->Leaves[Index_Leaf_y];
        sum_Delta_upto_Root(Table->Treeroot, Leaf, 
                            Pa->Event_Delta_Matrix[Type_of_Event+1][k]);
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Index_Node_y = Table->TOTAL_No_of_EVENTS * y + k;

        lambda_old   = Pa->rToI[k]; 
        time_old     = Table->Tree_Node_Index[Index_Node_y]->value;

        if(time_old == INFINITY) 
          assert(lambda_old == 0.0); 

        time_current = Rate->Stochastic_Time; 
      #endif

      Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event+1][k];

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Node = Table->Tree_Node_Index[Index_Node_y];
        lambda_new   = Pa->rToI[k];
        if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
          Pa->rToI[k] = 0.0;
          lambda_new  = 0.0; 
        }

        if (lambda_old == 0.0  && lambda_new == 0.0) 
          Next_Time = INFINITY;
        else if (lambda_old > 0.0   && lambda_new == 0.0)
          Next_Time = INFINITY;
        else if (lambda_old == 0.0  && lambda_new > 0.0)
          Next_Time = time_current - 1.0/Pa->rToI[k] * log(RANDOM);
        else if (lambda_old > 0.0  && lambda_new > 0.0) {
          if( time_old <= time_current ) {
            printf("time_old = %g\t time_current = %g\n", time_old, time_current);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) in Patches (%d, %d)\n", 
                    time_current, k, Type_of_Event, Patch[0], Patch[1]); 
            Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
            printtree(Table->Treeroot);

            assert( time_old > time_current );
          }
          Next_Time = time_current + lambda_old/lambda_new * (time_old-time_current);
        }
        else { 
          if(lambda_new < 0.0 || lambda_old < 0.0) {
            printf("A rate has become too negative: lambda_new = %g\t lambda_old = %g\n", 
                    lambda_new, lambda_old);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) in Patches (%d, %d)\n", 
                     time_current, k, Type_of_Event, Patch[0], Patch[1]);
            Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
            printtree(Table->Treeroot);
            Print_Press_Key(1,1,"Kill the program is lambda_new is too negative");
          }
        }      
        
        Node->value = Next_Time; 
        bubbling(Node, Table->Tree_Node_Index);
      #endif  
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
    // if( Rate->Total_Rate < -TOLERANCE) exit(0);
  }
}

void Updating_Event_Delta_Matrix(Community * Pa, int Type_of_Event, Parameter_Table * Table)
{
  double dC, Ga, ep; 
  double ** Delta_Matrix = Pa->Event_Delta_Matrix;
  int * n = Pa->n;
  double K_R = Table->K_R;
  double n_0, n_1; 
  double N_0 = (double)Pa->m_0;    /* Empty Space */

  n_0 = (double)n[0];
  n_1 = (double)n[1];

  assert( N_0 == (K_R - (n_0 + n_1)) );

  Ga = Table->Lambda_R_1 * Table->Chi_C_0;
  dC = 0.5 * Table->Lambda_R_1 * Table->p_2; 
  ep = Table->p_1;                      /* Segregation Error */ 

  switch(Type_of_Event) {
    case 0:
    case 2:
    case 3:
      /* The effect on system configuration of events (0), (2), and (3) is the same: 
                            n_0 ---> n_0 - 1
                            N_0 ---> N_0 + 1 
      Therefore, the change on the different rates is the same */
      Delta_Matrix[Type_of_Event][3] =  -dC * n_1/K_R;      
      Delta_Matrix[Type_of_Event][4] =  +Table->Beta_AP[0] * (n_0 - N_0 + 1.0)/K_R;
      Delta_Matrix[Type_of_Event][8] =  -dC * n_1/K_R;
      Delta_Matrix[Type_of_Event][9] =  +Table->Beta_AP[1] * (1.0 - ep) * n_1/K_R;
      Delta_Matrix[Type_of_Event][10] = +Table->Beta_AP[1] * ep * n_1/K_R;
      Delta_Matrix[Type_of_Event][11] = -Ga * n_1/K_R;  
      break;

    case 1:
    case 4:
    case 10:
      /* The effect on system configuration of events (1), (4), and (10) is the same: 
                            n_0 ---> n_0 + 1
                            N_0 ---> N_0 - 1 
      Therefore, the change on the different rates is the same */
      Delta_Matrix[Type_of_Event][3] =  +dC * n_1/K_R;      
      Delta_Matrix[Type_of_Event][4] =  +Table->Beta_AP[0] * (N_0 - n_0 + 1.0)/K_R;
      Delta_Matrix[Type_of_Event][8] =  +dC * n_1/K_R;
      Delta_Matrix[Type_of_Event][9] =  -Table->Beta_AP[1] * (1.0 - ep) * n_1/K_R;
      Delta_Matrix[Type_of_Event][10] = -Table->Beta_AP[1] * ep * n_1/K_R;
      Delta_Matrix[Type_of_Event][11] = +Ga * n_1/K_R;  
      break;

    case 5:
    case 7:
    case 8:
      /* The effect on system configuration of events (5), (7), and (8) is the same: 
                            n_1 ---> n_1 - 1
                            N_0 ---> N_0 + 1 
      Therefore, the change on the different rates is the same */
      Delta_Matrix[Type_of_Event][3] =  -dC * n_0/K_R;      
      Delta_Matrix[Type_of_Event][4] =  +Table->Beta_AP[0] * n_0/K_R;
      Delta_Matrix[Type_of_Event][8] =  -dC * n_0/K_R;
      Delta_Matrix[Type_of_Event][9] =  +Table->Beta_AP[1] * (1.0 - ep) * (n_1 - N_0 + 1.0)/K_R;
      Delta_Matrix[Type_of_Event][10] = +Table->Beta_AP[1] * ep * (n_1 - N_0 + 1.0)/K_R;
      Delta_Matrix[Type_of_Event][11] = -Ga * n_0/K_R;  
      break;

    case 6:
    case 9:
      /* The effect on system configuration of events (6) and (9) is the same: 
                            n_1 ---> n_1 + 1
                            N_0 ---> N_0 - 1 
      Therefore, the change on the different rates is the same */
      Delta_Matrix[Type_of_Event][3] =  +dC * n_0/K_R;      
      Delta_Matrix[Type_of_Event][4] =  -Table->Beta_AP[0] * n_0/K_R;
      Delta_Matrix[Type_of_Event][8] =  +dC * n_0/K_R;
      Delta_Matrix[Type_of_Event][9] =  +Table->Beta_AP[1] * (1.0 - ep) * (N_0 - n_1 + 1.0)/K_R;
      Delta_Matrix[Type_of_Event][10] = +Table->Beta_AP[1] * ep * (N_0 - n_1 + 1.0)/K_R;
      Delta_Matrix[Type_of_Event][11] = +Ga * n_0/K_R;  
      break;

    case 11:
      /* The effect on system configuration associated to a conjugation event (11) is: 
                            n_0 ---> n_0 - 1
                            n_1 ---> n_1 + 1
      Therefore, the change on the different rates is the same */
      Delta_Matrix[Type_of_Event][3] =  +dC * (n_0 - n_1 + 1.0 )/K_R;      
      Delta_Matrix[Type_of_Event][4] =  -Table->Beta_AP[0] * N_0/K_R;
      Delta_Matrix[Type_of_Event][8] =  +dC * (n_0 - n_1 + 1.0 )/K_R;
      Delta_Matrix[Type_of_Event][9] =  +Table->Beta_AP[1] * (1.0 - ep) * N_0/K_R;
      Delta_Matrix[Type_of_Event][10] = +Table->Beta_AP[1] * ep * N_0/K_R;
      Delta_Matrix[Type_of_Event][11] = +Ga * (n_0 - n_1 + 1.0 )/K_R;
      break;

    default:
      printf(" Error in Updating_Event_Delta_Matrix!!!\n");
      printf(" Type_of_Event = %d\t This value is not possible!!!\n", Type_of_Event);
      printf(" Only 0 to 11 are possible. Type of Event is ill-defined\n");
      printf(" The program will exit\n");
      Print_Press_Key(1,0,"."); 
      exit(0);
  }
}
