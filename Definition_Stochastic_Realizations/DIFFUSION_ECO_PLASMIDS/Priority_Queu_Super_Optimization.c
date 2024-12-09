#include <MODEL.h>

#define TOLERANCE 1.0E-9

/* This functions allocate, initialize and free a number of local communities,
   which make up our total patch system or metapopulation */
extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Priority_Queu_Super_Optimization( Community * Pa, int x, Parameter_Table * Table, 
                                       int n, int Sp, int k , double time_current,  
                                       double lambda_old, double lambda_new )
{
  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double time_old, Next_Time;
  treenode * Leaf;
  treenode * Node;

  Index_Node_x = Table->TOTAL_No_of_EVENTS * x + n;

  time_old     = Table->Tree_Node_Index[Index_Node_x]->value;
                
  if(time_old == INFINITY) assert(lambda_old == 0.0);
    
  Node = Table->Tree_Node_Index[Index_Node_x];
    
          if (lambda_new > -TOLERANCE && lambda_new < 0.0 ) { 
            Pa->rToI[n] = 0.0;
            lambda_new  = 0.0; 
          }

          if(k == Type_Event) {
            if( Pa->rToI[n] == 0.0 ) 
              Next_Time = INFINITY;
            else
              Next_Time = time_current - 1.0/Pa->rToI[n] * log(RANDOM);
                                  
            ( * bool_Next_Time) = true; 
          }
          else if (lambda_old == 0.0  && lambda_new == 0.0) 
            Next_Time = INFINITY;
          else if (lambda_old > 0.0   && lambda_new == 0.0)
            Next_Time = INFINITY;
          else if (lambda_old == 0.0  && lambda_new > 0.0)
            Next_Time = time_current - 1.0/Pa->rToI[n] * log(RANDOM);
          else if (lambda_old > 0.0  && lambda_new > 0.0) {
            if( time_old <= time_current ) {
              printf("time_old = %g\t time_current = %g\n", time_old, time_current);
              printf("Time = %g\t Related Event, from 0 to %d (or conjugating species pair Sp_0) = %d(%d) affecting Sp %d in Patches (%d, %d)\n", 
                      time_current, Table->No_of_EVENTS-1, k, Type_Event, Sp, Patch[0], Patch[1]); 
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
} 

