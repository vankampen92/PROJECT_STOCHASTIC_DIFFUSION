#include <MODEL.h>

#define TOLERANCE 1.0E-9

extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Temporal_Dynamics_Update_Focal_Species ( Community * Pa, int x, Parameter_Table * Table, 
                                              Stochastic_Rate * Rate, 
                                              int Type_of_Event, int Type_Event, 
                                              bool * bool_Next_Time, double * delta_rate )
{
  /* Input: 
      . Pa, focal patch 
           
     Output: 
      . Delta_Rate, 
  */ 
  int i, k, x, n, m, Sp, Sp_0;

  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double time_old, lambda_old, lambda_new, Next_Time;
  treenode * Leaf;
  treenode * Node;  

  double time_current = Rate->Stochastic_Time;
  double Delta_Rate = * delta_rate; 

  Sp         = Type_of_Event/(Table->No_of_EVENTS-1);
  
  Updating_Event_Delta_Matrix_Partial_Events(Pa, Type_Event, Sp, Table); 
        
  m = Pa->Event_Adjacence_List[Type_Event][n];   /* How many events are connected 
                                                    to event 'Type_Event'???
                                                    i.e., the lenght of the adjacence  
                                                    list of 'Type_of_Event'
                                                 */
  for(i=0; i < (m-1); i++) {
    k = Pa->Event_Adjacence_List[Type_Event][i];  /* Which events are connected 
                                                     to event 'Type_of_Event'???
                                                  */
    Delta_Rate  += Pa->Event_Delta_Tensor[Sp][Type_Event][k];

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
          Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Sp*(Table->No_of_EVENTS-1) + k;
          Leaf = Table->Leaves[Index_Leaf_x];
          sum_Delta_upto_Root(Table->Treeroot, Leaf, Pa->Event_Delta_Tensor[Sp][Type_Event][k]);
    #endif                            

    lambda_old                                = Pa->rToI[Sp*(Table->No_of_EVENTS-1) + k];  
    Pa->rToI[Sp*(Table->No_of_EVENTS-1) + k] += Pa->Event_Delta_Tensor[Sp][Type_Event][k];
    lambda_new                                = Pa->rToI[Sp*(Table->No_of_EVENTS-1) + k];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        n = Sp*(Table->No_of_EVENTS-1) + k;  
        Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp, k, time_current, 
                                          lambda_old, lambda_new );                   
    #endif
  }

  #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      if( (* bool_Next_Time) == false) {
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
  
  * delta_rate = Delta_Rate; 
}

