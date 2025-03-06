/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2025 (c)                        */
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
                    Patch[2] is the species involved in the event.
                    Patch[3] is the other species involved in the event.
                    Patch[4] is a flag
     Output arguments: 
     . Rate         Stochastic Rate is updated from previous value (without recalculating)
  */
  int i, n, m, k, Sp; 
  int x, y; 
  
  Community * Pa;

  n = Table->TOTAL_No_of_EVENTS;
  
  x = Patch[0]; y = Patch[1]; Sp = Patch[2]; /* First Species involved */

  if (Type_of_Event >= n || Type_of_Event < 0 ) {
    printf(" Error in Temporal Dynamics Update!!!\n");
    printf(" Type of Event occurring is too large.\n");
    printf(" It can only be labeled from 0 or %d\n", Table->TOTAL_No_of_EVENTS-1);
    printf(" but events is seemingly labeled as %d\n", Type_of_Event); 
    printf(" The program will exit\n"); 
    exit(0); 
  }

  if ( x == y ) {  /* Local event on a single Patch x */
  
    assert( Type_of_Event != 0 && Type_of_Event != 5 ) ; /* Because (0) and (5) are 
                                                            out migration events. Only
							    possible when x patch 
							    is different from y patch
							 */
    Pa    = My_Community[x];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event, Table); 
    Update_Event_Rate_Structure( Pa, x, Type_of_Event, Sp, Rate, Table, 0 );

  }
  else {  /* MOVEMENT EVENT: x: patch exporting an individual  
                             y: patch receiving an individual, 
             0: Plasmid Free Cells at rate Table->Mu
             1: Plasmid Carrying Cells at rate Table->Mu_C */  		      
    assert( Type_of_Event == 0 || Type_of_Event == 5 ) ; 
    assert( Table->No_of_CELLS > 1 );
    
    /* Changes in rates due to the loss of an individual in patch x */
    Pa    = My_Community[x];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event, Table);
    Update_Event_Rate_Structure( Pa, x, Type_of_Event, Sp, Rate, Table, 0 ); 

    /* Changes in rates due to the adquisition of an individual in patch y */
    Pa    = My_Community[y];
    Updating_Event_Delta_Matrix(Pa, Type_of_Event+1, Table);
    Update_Event_Rate_Structure( Pa, y, Type_of_Event+1, Sp, Rate, Table, 1 );
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

void Update_Event_Rate_Structure( Community * Pa, int x, int Type_of_Event, int Sp,
				                          Stochastic_Rate * Rate, Parameter_Table * Table,
				                          int patch)
{
  bool bool_Next_Time;
  int i, k, m, n;
  double Delta_Rate;
  double time_old, time_current, Next_Time; 
  double lambda_old, lambda_new;
  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  treenode * Leaf;
  treenode * Node;

  /* Input arguments: 
     . x, Integer label of the local patch where the event has occurred. 
     . Pa, Pointer to the patch itself. 
     . Sp, Integer label of the species affected. 
     . Rate, Stochastic_Rate containing the current time (time_current)
     . Table, Parameter Table Structure, as usual
       . patch should be zero most of the time. Only patch = 1 in 
       in a movement event when updating the rates at the 2nd patch
  */

  // #if defined BINARY_TREE_SUPER_OPTIMIZATION
    /* The index of the tip the tree corresponding to a particular 
       event in a given cell 
    */
  //   Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Type_of_Event;  
  // #endif 

  n = Table->TOTAL_No_of_EVENTS; 
  m = Pa->Event_Adjacence_List[Type_of_Event][n]; /* How many events are affected by 
						                                         the configuration change induced
						                                         by the event 'Type_of_Event'???
						                                         i.e., the lenght of the adjacence  
						                                         list of 'Type_of_Event'
						                                      */
  time_current = Rate->Stochastic_Time;
  bool_Next_Time = false; 
  Delta_Rate = 0.0;
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
    time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
    lambda_old   = Pa->rToI[k]; 
    
    if(time_old == INFINITY) {
      if(lambda_old != 0.0) {
	      printf("time_old = %g\t time_current = %g\n", 
	              time_old, time_current);
	      printf("lambda_old = %g\t lambda_new = %g\n", 
	              lambda_old, Pa->rToI[k]);
	      printf("Related Event (Type of Event) = %d(%d) in Patch  %d\n",
	              k, Type_of_Event, x);        
	      assert(lambda_old == 0.0);
      }
    }
#endif 
    
    Pa->rToI[k] += Pa->Event_Delta_Matrix[Type_of_Event][k];
    
#if defined PRIORITY_QUEU_SUPER_OPTIMIZATION 
    Priority_Queu_Super_Optimization( Pa, x, Type_of_Event, Sp, k,
				                              time_current, lambda_old,
				                              Table,
				                              &bool_Next_Time);
    // Print_Press_Key(1,1,"Printing out Tree before after bubbling\n");
    // printtree(Table->Treeroot);
#endif  
  }

#if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
    if( bool_Next_Time == false && patch == 0) {
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
    
void Priority_Queu_Super_Optimization( Community * Pa, int x, 
                                       int n, int Sp, int k ,
				                               double time_current, 
                                       double lambda_old,
				                               Parameter_Table * Table,
				                               bool * bool_Next_Time)
{
  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double time_old, Next_Time;
  double lambda_new;
  treenode * Leaf;
  treenode * Node;
  /* This function performs the updating of the Priority Queu stored in the 
     Binary Tree as a consequence of the event that has just occurred (n), and, therefore, 
     acording to the configurational change generated. Only the the (k) event rate is
     updated (see input argument k), its next time is re-calculated and bubbled finally   
     in the bubbling algorithm to keep the priority order up to date. The event (k) should 
     be one the the events whose rate is being affected by the configurational change 
     induced by the event (n). The rest of affected events are stored in the coreesponding
     row of the adjancence list of the event (n) (i.e., Pa->Event_Adjacence_List[n][] array).
     
     Input parameters: 

     . n, Integer label of the event occurred (Type_of_Event).
     . x, Integer label of the local patch where the event has occurred. 
     . Pa, Pointer to the patch itself. 
     . Sp, Integer label of the species affected. 
     . k,  Integer label of an event whose rate is affected by the 
           configurational change induced by the event 'n' that has 
           just occurred. 
     . time_current 
     . lambda_old,  value of old rate of event k 
     . lambda_new,  value (pointer to) of the new rate of event k  
     . Table, Parameter Table Structure, as usual
  */
  
  Index_Node_x = Table->TOTAL_No_of_EVENTS * x + k;

  time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
               
  Node         = Table->Tree_Node_Index[Index_Node_x];

  if(time_old == INFINITY) assert(lambda_old == 0.0);

  lambda_new   = Pa->rToI[k];

  if (lambda_new < 0.0 ) { 
    if(lambda_new < -TOLERANCE ) {
      printf("Warning: Rates might be too negative. This may indicate an error!!!\n");
      printf("lambda_new = %g\t lambda_old = %g\n", lambda_new, lambda_old);
      printf("Time = %g\t Event (between 0 to %d) = %d, whose rate has been\n", 
              time_current, Table->TOTAL_No_of_EVENTS-1, k);
      printf("affected by the execution of event %d (affecting Sp %d) in Patch %d\n",  
	            n, Sp, x);  
    }

    Pa->rToI[k] = 0.0;
    lambda_new  = 0.0;  
  }

  if(k == n) {
    if( Pa->rToI[k] == 0.0 )
      Next_Time = INFINITY;
    else {
       Next_Time = time_current - 1.0/Pa->rToI[k] * log(RANDOM);                        
       ( * bool_Next_Time) = true;
    }
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
      printf("Time = %g\t Event (between 0 to %d) = %d, whose rate has been\n", 
              time_current, Table->TOTAL_No_of_EVENTS-1, k);
      printf("affected by the execution of event %d (affecting Sp %d) in Patch %d\n",  
	            n, Sp, x); 
      Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");
      
      if(Table->TOTAL_GRAND_No_of_EVENTS <= 16) printtree(Table->Treeroot);
      assert( time_old > time_current );
    }
    Next_Time = time_current + lambda_old/lambda_new * (time_old-time_current);
  }
  else { 
    if(lambda_new < 0.0 || lambda_old < 0.0) {
      printf("A rate has become too negative: lambda_new = %g\t lambda_old = %g\n", 
	            lambda_new, lambda_old);
      printf("Time = %g\t Event (between 0 to %d) = %d, whose rate has been\n", 
              time_current, Table->TOTAL_No_of_EVENTS-1, k);
      printf("affected by the execution of event %d (affecting Sp %d) in Patch %d\n",  
              n, Sp, x); 
      Print_Press_Key(1,1,"Printing out Tree in Temporal_Dynamics_Update...");

      if(Table->TOTAL_GRAND_No_of_EVENTS <= 16) printtree(Table->Treeroot);
      
      Print_Press_Key(1,1,"Kill the program is lambda_new is too negative");
      
      if ( MIN(lambda_new, lambda_old)  < (-TOLERANCE) ) {
	      printf("Rates too negative. The program will exit\n");
	      exit(0); 
      }
    }
  }

  Node->value = Next_Time; 
  bubbling(Node, Table->Tree_Node_Index);
  // Print_Press_Key(1,1,"Printing out Tree before after bubbling\n");
  // printtree(Table->Treeroot);
} 
