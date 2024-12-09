/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2021 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

#define TOLERANCE 1.0E-9

void Updating_Event_Delta_Matrix(Community * Pa, int Type_of_Event, int Sp, int Sp_Out,  Parameter_Table * Table);

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
     . Type_of_Event is an index to the local event occurring in a Patch between 0 and Table->TOTAL_No_of_EVENTS-1
     . Patch        is an array containing the two patches involved in a movement event. 
                    Patch[0] is the patch sending the individual
		                Patch[1] is the patch receiving the individual.
                    Patch[2] is the index of the species receiving a mutant

     Output arguments: 
     . Rate         Stochastic Rate is updated from previous value (without recalculating)
  */
  int Type_Event;      /* Index of the event between 0 and Table->No_of_EVENTS-1 */
                       /* Type of Event that occurs at a given species           */
  bool bool_Next_Time;
  int i, n, m, k; 
  int x, y, Sp, Sp_Out; 
  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double Delta_Rate, time_old, lambda_old, lambda_new, time_current, Next_Time;
  treenode * Leaf;
  treenode * Node;  
  
  Community * Pa;
  double OutMigration; 

  n = Table->No_of_EVENTS;
  
  x = Patch[0]; y = Patch[1]; Sp_Out = Patch[2];

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

    Sp         = Type_of_Event/Table->No_of_EVENTS;
    Type_Event = Type_of_Event%Table->No_of_EVENTS;
  
    assert( Type_Event != 0 ) ; /* Because this is an out
							                     migration events. Only
							                     possible when x patch 
							                     is different from y patch
							                  */
    Pa    = My_Community[x];
       
    Updating_Event_Delta_Matrix(Pa, Type_Event, Sp, Sp_Out, Table); 
        
    m = Pa->Event_Adjacence_List[Type_Event][n];   /* How many events are connected 
                                                        to event 'Type_Event'???
                                                        i.e., the lenght of the adjacence  
                                                        list of 'Type_of_Event'
                                                     */
    Delta_Rate = 0.0;
    bool_Next_Time = false; 
    time_current = Rate->Stochastic_Time; 
    for(i=0; i<m; i++) {
        k = Pa->Event_Adjacence_List[Type_Event][i]; /* Which events are connected 
                                                          to event 'Type_of_Event'???
                                                        */
        Delta_Rate  += Pa->Event_Delta_Tensor[Sp][Type_Event][k];

        #if defined BINARY_TREE_SUPER_OPTIMIZATION
          Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Sp*Table->No_of_EVENTS + k;
          Leaf = Table->Leaves[Index_Leaf_x];
          sum_Delta_upto_Root(Table->Treeroot, Leaf, Pa->Event_Delta_Tensor[Sp][Type_Event][k]);
        #endif                            

        #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
          Index_Node_x = Table->TOTAL_No_of_EVENTS * x + Sp*Table->No_of_EVENTS + k;

          lambda_old   = Pa->rToI[Sp*Table->No_of_EVENTS + k]; 
          time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
                
          if(time_old == INFINITY) assert(lambda_old == 0.0);
        #endif 

        Pa->rToI[Sp*Table->No_of_EVENTS + k] += Pa->Event_Delta_Tensor[Sp][Type_Event][k];

        #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
          Node = Table->Tree_Node_Index[Index_Node_x];
          lambda_new   = Pa->rToI[Sp*Table->No_of_EVENTS + k];
          if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
            Pa->rToI[Sp*Table->No_of_EVENTS + k] = 0.0;
            lambda_new  = 0.0; 
          }

          if(k == Type_Event) {
            if( Pa->rToI[Sp*Table->No_of_EVENTS + k] == 0.0 ) 
              Next_Time = INFINITY;
            else
              Next_Time = time_current - 1.0/Pa->rToI[Sp*Table->No_of_EVENTS + k] * log(RANDOM);
                                  
            bool_Next_Time = true; 
          }
          else if (lambda_old == 0.0  && lambda_new == 0.0) 
            Next_Time = INFINITY;
          else if (lambda_old > 0.0   && lambda_new == 0.0)
            Next_Time = INFINITY;
          else if (lambda_old == 0.0  && lambda_new > 0.0)
            Next_Time = time_current - 1.0/Pa->rToI[Sp*Table->No_of_EVENTS + k] * log(RANDOM);
          else if (lambda_old > 0.0  && lambda_new > 0.0) {
            if( time_old <= time_current ) {
              printf("time_old = %g\t time_current = %g\n", time_old, time_current);
              printf("Time = %g\t Related Event (Type Event) = %d(%d) affecting Sp %d in Patches (%d, %d)\n", 
                      time_current, k, Type_Event, Sp, Patch[0], Patch[1]); 
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
    /* Extra affectation of rates from other species/types when Event is the
       establishment of an adult plant (4 or 5) */
    if (Type_Event == 5 ) {
      /* Changes in rates due to the establisment of a mutant as adult platn */ 
      assert( Sp_Out < Table->No_of_RESOURCES); 

      Updating_Event_Delta_Matrix(Pa, Type_Event+1, Sp, Sp_Out, Table);
    
      m = Pa->Event_Adjacence_List[Type_Event+1][n]; /* How many events are connected 
							                                          to event 'Type_Event+1'???
							                                          m is the lenght of the adjacence list 
							                                          of 'Type_Event+1', 
							                                          because 6 is an event  
							                                          involving the increase of 
							                                          the plant established population 
                                                        of the Sp_Out mutated species 
						                                        */
      for(i=0; i<m; i++) {
        k = Pa->Event_Adjacence_List[Type_Event+1][i]; /* Which events are connected 
                                                          to event 'Type_of_Event'???
                                                        */
        Delta_Rate  += Pa->Event_Delta_Tensor[Sp_Out][Type_Event+1][k];

        #if defined BINARY_TREE_SUPER_OPTIMIZATION
          Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Sp_Out*Table->No_of_EVENTS + k;
          Leaf = Table->Leaves[Index_Leaf_x];
          sum_Delta_upto_Root(Table->Treeroot, Leaf, Pa->Event_Delta_Tensor[Sp_Out][Type_Event+1][k]);
        #endif                            

        #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
          Index_Node_x = Table->TOTAL_No_of_EVENTS * x + Sp_Out*Table->No_of_EVENTS + k;

          lambda_old   = Pa->rToI[Sp_Out*Table->No_of_EVENTS + k]; 
          time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
                
          if(time_old == INFINITY) assert(lambda_old == 0.0);
        #endif 

        Pa->rToI[Sp_Out*Table->No_of_EVENTS + k] += Pa->Event_Delta_Tensor[Sp_Out][Type_Event+1][k];

        #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
          Node = Table->Tree_Node_Index[Index_Node_x];
          lambda_new   = Pa->rToI[Sp_Out*Table->No_of_EVENTS + k];
          if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
            Pa->rToI[Sp_Out*Table->No_of_EVENTS + k] = 0.0;
            lambda_new  = 0.0; 
          }

          if (lambda_old == 0.0  && lambda_new == 0.0) 
            Next_Time = INFINITY;
          else if (lambda_old > 0.0   && lambda_new == 0.0)
            Next_Time = INFINITY;
          else if (lambda_old == 0.0  && lambda_new > 0.0)
            Next_Time = time_current - 1.0/Pa->rToI[Sp_Out*Table->No_of_EVENTS + k] * log(RANDOM);
          else if (lambda_old > 0.0  && lambda_new > 0.0) {
            if( time_old <= time_current ) {
              printf("time_old = %g\t time_current = %g\n", time_old, time_current);
              printf("Time = %g\t Related Event (Type Event) = %d(%d) affecting Sp %d in Patches (%d, %d)\n", 
                      time_current, k, Type_Event, Sp_Out, Patch[0], Patch[1]); 
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
                      time_current, k, Type_Event, Patch[0], Patch[1]);
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
    }
    if( Type_Event == 4 || Type_Event == 5 || Type_Event == 7) {
      /* Establisment of a new adult plant with mutatation (Type_Event = 5) 
         or without mutation (Type_Event = 4), 
         or death of an adult plant (which opens new space, Type_Event = 7) 
         In any of these cases, the rates of establishment of the other species also 
         change 
      */
      for(i=0; i<Table->No_of_RESOURCES; i++) {
        if( i != Sp && i != Sp_Out )  {
          for (k = 4; k<=5; k++) {
            /* k will only take 2 values (4, 5) */ 
            Delta_Rate  += Pa->Event_Delta_Tensor[i][Type_Event][k];
            #if defined BINARY_TREE_SUPER_OPTIMIZATION
              Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + i*Table->No_of_EVENTS + k;
              Leaf = Table->Leaves[Index_Leaf_x];
              sum_Delta_upto_Root(Table->Treeroot, Leaf, Pa->Event_Delta_Tensor[i][Type_Event][k]);
            #endif                            

            #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
              Index_Node_x = Table->TOTAL_No_of_EVENTS * x + i*Table->No_of_EVENTS + k;

              lambda_old   = Pa->rToI[i*Table->No_of_EVENTS + k]; 
              time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
                
              if(time_old == INFINITY) assert(lambda_old == 0.0);
            #endif 

            Pa->rToI[i*Table->No_of_EVENTS + k] += Pa->Event_Delta_Tensor[i][Type_Event][k];

            #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
              Node = Table->Tree_Node_Index[Index_Node_x];
              lambda_new   = Pa->rToI[i*Table->No_of_EVENTS + k];
              if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
                Pa->rToI[i*Table->No_of_EVENTS + k] = 0.0;
                lambda_new  = 0.0; 
              }

              if (lambda_old == 0.0  && lambda_new == 0.0) 
                Next_Time = INFINITY;
              else if (lambda_old > 0.0   && lambda_new == 0.0)
                Next_Time = INFINITY;
              else if (lambda_old == 0.0  && lambda_new > 0.0)
                Next_Time = time_current - 1.0/Pa->rToI[i*Table->No_of_EVENTS + k] * log(RANDOM);
              else if (lambda_old > 0.0  && lambda_new > 0.0) {
                if( time_old <= time_current ) {
                  printf("time_old = %g\t time_current = %g\n", time_old, time_current);
                  printf("Time = %g\t Related Event (Type Event) = %d(%d) affecting Sp %d in Patches (%d, %d)\n", 
                      time_current, k, Type_Event, i, Patch[0], Patch[1]); 
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
        }
      }
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
  else {  /* MOVEMENT EVENT involving two Patches: 
             x: patch exporting an individual propagule 
             y: patch receiving an individual propagule, 
             out from patch 'x' into patch 'y'   
          */       
		      /* Out migration event sending one individual:
             RP: propagules
          */
    Sp         = Type_of_Event/Table->No_of_EVENTS;
    Type_Event = Type_of_Event%Table->No_of_EVENTS;         

    assert( Type_Event == 0 ) ; 
    assert( Table->No_of_CELLS > 1 );
    
    /* Changes in rates in patch x due to out movement 
       of one individual from patch x to patch y  
    */
    Pa    = My_Community[x];
    Updating_Event_Delta_Matrix(Pa, Type_Event, Sp, Sp_Out, Table);
    
    m = Pa->Event_Adjacence_List[Type_Event][n]; /* How many events are connected 
						                                        to event 'Type_of_Event'???
						                                        m is the lenght of the adjacence list 
						                                        of 'Type_of_Event' */
    Delta_Rate = 0.0;
    bool_Next_Time = false;
    time_current = Rate->Stochastic_Time;
    for(i=0; i<m; i++) {
      k = Pa->Event_Adjacence_List[Type_Event][i]; /* Which events are connected 
							                                        to event 'Type_of_Event'???
						                                        */
      Delta_Rate  += Pa->Event_Delta_Tensor[Sp][Type_Event][k];
      
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Sp*Table->No_of_EVENTS + k;
        Leaf = Table->Leaves[Index_Leaf_x];
        sum_Delta_upto_Root(Table->Treeroot, Leaf, Pa->Event_Delta_Tensor[Sp][Type_Event][k]);
      #endif

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Index_Node_x = Table->TOTAL_No_of_EVENTS * x + Sp*Table->No_of_EVENTS + k;

        lambda_old   = Pa->rToI[Sp*Table->No_of_EVENTS + k]; 
        time_old     = Table->Tree_Node_Index[Index_Node_x]->value;

        if(time_old == INFINITY) assert(lambda_old == 0.0);
      #endif 

      Pa->rToI[Sp*Table->No_of_EVENTS + k] += Pa->Event_Delta_Tensor[Sp][Type_Event][k];

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Node = Table->Tree_Node_Index[Index_Node_x];
        lambda_new   = Pa->rToI[k];
        if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
          Pa->rToI[Sp*Table->No_of_EVENTS + k] = 0.0;
          lambda_new  = 0.0; 
        }

        if(k == Type_Event) {
          if( Pa->rToI[Sp*Table->No_of_EVENTS + k] == 0.0 ) 
            Next_Time = INFINITY;
          else
            Next_Time = time_current - 1.0/Pa->rToI[Sp*Table->No_of_EVENTS + k] * log(RANDOM);
                                
          bool_Next_Time = true; 
        }
        else if (lambda_old == 0.0  && lambda_new == 0.0) 
          Next_Time = INFINITY;
        else if (lambda_old > 0.0   && lambda_new == 0.0)
          Next_Time = INFINITY;
        else if (lambda_old == 0.0  && lambda_new > 0.0)
          Next_Time = time_current - 1.0/Pa->rToI[Sp*Table->No_of_EVENTS + k] * log(RANDOM);
        else if (lambda_old > 0.0  && lambda_new > 0.0) {
          if( time_old <= time_current ) {
            printf("time_old = %g\t time_current = %g\n", time_old, time_current);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) affecting Sp %d in Patches (%d, %d)\n", 
                    time_current, k, Type_of_Event, Sp, Patch[0], Patch[1]); 
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
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) affecting Sp %d in Patches (%d, %d)\n", 
                     time_current, k, Type_of_Event, Sp, Patch[0], Patch[1]);
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

    /* Changes in rates due to 
       the adquisition of an individual propagule of species Sp in patch y */
    Pa    = My_Community[y];
    Updating_Event_Delta_Matrix(Pa, Type_Event+1, Sp, Sp_Out, Table);
    
    m = Pa->Event_Adjacence_List[Type_Event+1][n]; /* How many events are connected 
							                                        to event 'Type_Event+1'???
							                                        m is the lenght of the adjacence list 
							                                        of 'Type_Event+1', 
							                                        because 1 is an event  
							                                        involving the increase of 
							                                        the propagule population 
						                                        */
    Delta_Rate = 0.0;
    for(i=0; i<m; i++) {
      k = Pa->Event_Adjacence_List[Type_Event+1][i]; /* Which events are connected 
							                                             to event 'Type_of_Event'???
							                                          */
      Delta_Rate  += Pa->Event_Delta_Tensor[Sp][Type_Event+1][k];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Index_Leaf_y = Table->TOTAL_No_of_EVENTS * y + Sp*Table->No_of_EVENTS + k;  
        Leaf = Table->Leaves[Index_Leaf_y];
        sum_Delta_upto_Root(Table->Treeroot, Leaf, 
                            Pa->Event_Delta_Tensor[Sp][Type_Event+1][k]);
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Index_Node_y = Table->TOTAL_No_of_EVENTS * y + Sp*Table->No_of_EVENTS + k;

        lambda_old   = Pa->rToI[Sp*Table->No_of_EVENTS + k]; 
        time_old     = Table->Tree_Node_Index[Index_Node_y]->value;

        if(time_old == INFINITY) 
          assert(lambda_old == 0.0); 

        time_current = Rate->Stochastic_Time; 
      #endif

      Pa->rToI[Sp*Table->No_of_EVENTS + k] += Pa->Event_Delta_Tensor[Sp][Type_Event+1][k];

      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Node = Table->Tree_Node_Index[Index_Node_y];
        lambda_new   = Pa->rToI[Sp*Table->No_of_EVENTS + k];
        if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
          Pa->rToI[Sp*Table->No_of_EVENTS + k] = 0.0;
          lambda_new  = 0.0; 
        }

        if (lambda_old == 0.0  && lambda_new == 0.0) 
          Next_Time = INFINITY;
        else if (lambda_old > 0.0   && lambda_new == 0.0)
          Next_Time = INFINITY;
        else if (lambda_old == 0.0  && lambda_new > 0.0)
          Next_Time = time_current - 1.0/Pa->rToI[Sp*Table->No_of_EVENTS + k] * log(RANDOM);
        else if (lambda_old > 0.0  && lambda_new > 0.0) {
          if( time_old <= time_current ) {
            printf("time_old = %g\t time_current = %g\n", time_old, time_current);
            printf("Time = %g\t Related Event (Type of Event) = %d(%d) affecting Species %d in Patches (%d, %d)\n", 
                    time_current, k, Type_of_Event, Sp, Patch[0], Patch[1]); 
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

void Updating_Event_Delta_Matrix(Community * Pa, int Type_of_Event, int Sp, int Sp_Out, Parameter_Table * Table)
{
  /* 
     This is the subset of the Delta Matrix entries that depend
     on system configuration. Therefore, they need to be 
     changed in agreement to the event that has just occurred 
  */
  int i; 

  double ** Delta_Matrix = Pa->Event_Delta_Tensor[Sp];

  int * n = Pa->n;

  double K_R = Table->K_R;
  int m_0    = Pa->m_0;

  /* Definition of the local state vector numerical order */
  #include <Model_Variables_Code.Include.c>

   RP = 2*Sp + Table->RP; 
   R  = 2*Sp + Table->R; 
  
  switch( Type_of_Event )
    {
    case 0:  /* 0: Propagule Out-Migration (P --> P-1) and some other patch gains one */ 
      Delta_Matrix[0][4] = -Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (double)m_0; 
      Delta_Matrix[0][5] = -Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (double)m_0;
    break;
      
    case 1:  /* 1: Propagule External Immigration event */
      Delta_Matrix[1][4] = +Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (double)m_0; 
      Delta_Matrix[1][5] = +Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (double)m_0;
      
    break;
    
    case 2: /* 2: Propagule Death  */
      Delta_Matrix[2][4] = -Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (double)m_0; 
      Delta_Matrix[2][5] = -Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (double)m_0;
      
    break;
    
    case 3: /* 3: Propagule Production */
      Delta_Matrix[3][4] = +Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (double)m_0; 
      Delta_Matrix[3][5] = +Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (double)m_0;
      
    break;
      
    case 4: /* 4: Propagule Establishment  (withtout mutation)*/
      RP = 2*Sp + Table->RP; 
        
      Delta_Matrix[4][4] = -Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (1.0 + (double)n[RP] + (double)m_0);
      Delta_Matrix[4][5] = -Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (1.0 + (double)n[RP] + (double)m_0);

      for(i=0; i<Table->No_of_RESOURCES; i++)
        if( i != Sp) {
          Delta_Matrix = Pa->Event_Delta_Tensor[i];
          RP = 2*i + Table->RP; 
          
          Delta_Matrix[4][4] = -Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (double)n[RP];
          Delta_Matrix[4][5] = -Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (double)n[RP];
        }     
      
      Delta_Matrix = Pa->Event_Delta_Tensor[Sp];
    break;
    
    case 5:  /* 5: Propagule Establishment (with mutatation)  */
      RP = 2*Sp + Table->RP; 
        
      Delta_Matrix[5][4] = -Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (1.0 + (double)n[RP] + (double)m_0);
      Delta_Matrix[5][5] = -Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (1.0 + (double)n[RP] + (double)m_0);

      for(i=0; i<Table->No_of_RESOURCES; i++)
        if( i != Sp) {
          Delta_Matrix = Pa->Event_Delta_Tensor[i];
          RP = 2*i + Table->RP; 
          
          Delta_Matrix[5][4] = -Table->Eta_RP[i]*(1.0-2.0*Table->p_1)/K_R * (double)n[RP];
          Delta_Matrix[5][5] = -Table->Eta_RP[i]*(2.0*Table->p_1)    /K_R * (double)n[RP];
        }     
      
      Delta_Matrix = Pa->Event_Delta_Tensor[Sp];

    break;
    
    case 6:  /* 6: Mutation in (in Species Sp_Out) */ 
      Delta_Matrix = Pa->Event_Delta_Tensor[Sp_Out];
      RP = 2*Sp_Out + Table->RP; 
          
      Delta_Matrix[6][4] = -Table->Eta_RP[Sp_Out]*(1.0-2.0*Table->p_1)/K_R * (double)n[RP];
      Delta_Matrix[6][5] = -Table->Eta_RP[Sp_Out]*(2.0*Table->p_1)    /K_R * (double)n[RP];

      Delta_Matrix = Pa->Event_Delta_Tensor[Sp];
      
    break;

    case 7:  /* 7: Adult plant death (AP ---> AP - 1) */  
      Delta_Matrix = Pa->Event_Delta_Tensor[Sp];
      RP = 2*Sp + Table->RP; 
          
      Delta_Matrix[7][4] = +Table->Eta_RP[Sp]*(1.0-2.0*Table->p_1)/K_R * (double)n[RP];
      Delta_Matrix[7][5] = +Table->Eta_RP[Sp]*(2.0*Table->p_1)    /K_R * (double)n[RP];
      
      for(i=0; i<Table->No_of_RESOURCES; i++)
        if( i != Sp) {
          Delta_Matrix = Pa->Event_Delta_Tensor[i];
          RP = 2*i + Table->RP; 
          
          Delta_Matrix[7][4] = +Table->Eta_RP[i]*(1.0-2.0*Table->p_1)/K_R * (double)n[RP];
          Delta_Matrix[7][5] = +Table->Eta_RP[i]*(2.0*Table->p_1)    /K_R * (double)n[RP];
        }     
      
      Delta_Matrix = Pa->Event_Delta_Tensor[Sp];

    break;

    default:
      /* Something is very very wrong!!! */
      printf(" Type_of_Event = %d\t This value is not possible!!!\n", Type_of_Event);
      printf(" Only 0 to 7 are possible. Type of Event is ill-defined\n");
      printf(" The program will exit\n");
      Print_Press_Key(1,0,"."); 
      exit(0);
    }
}
