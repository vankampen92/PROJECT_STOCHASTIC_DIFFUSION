#include <MODEL.h>

#define TOLERANCE 1.0E-9

extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Temporal_Dynamics_Update_Competition_Rates ( Community * Pa, int x, Parameter_Table * Table, 
                                                  Stochastic_Rate * Rate, 
                                                  int Type_of_Event, int Type_Event, 
                                                  bool * bool_Next_Time, double * delta_rate )
{
    /* Changing rates of the competition-induced death of the species competing with
       the focal one (in ist competition list).  Updating Per Capita Competition rates: required!!!. 
    */ 
    /* Input: 
      . Pa, focal patch 
      . x
      . Table
      . Rate
      . Type_of_Event
      . Type_Event

      Output: 
      . Delta_Rate.
      . bool_Next_Time. 
    */ 
  double f; 
  int i, k, x, n, Sp, Sp_0;
  double Delta_Vector; 
  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double time_old, lambda_old, lambda_new, Next_Time;
  treenode * Leaf;
  treenode * Node;   

  double time_current = Rate->Stochastic_Time;
  double Delta_Rate = * delta_rate; 
       
  if(Type_Event == 0 || Type_Event == 2 || Type_Event == 3)
        f = (-1.0);
  else if (Type_Event == 1 || Type_Event == 4)
        f = ( 1.0);
  else {
        printf("Error in Temporal_Dynamics_Empty_Conjugation_Rates_Focal_Species.c\n");
        exit(0); 
  }

  Sp_0  = Type_of_Event/(Table->No_of_EVENTS-1); /* Focal Species */
  
  /* Conjugation Rates in Competition List of the Focal Species */
  for( j = 0; j < Table->Competition_List_Indeces[Sp_0][Table->No_of_RESOURCES]; j++) { 
    
    Sp = Table->Competition_List_Indeces[Sp_0][j]; /* Strain ID of a recipient of species Sp_0 */
    assert( Sp != Sp_0); 

    Delta_Vector = f * Table->Competition_Induced_Death[Sp_0][j] / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;

    k = 3; /* Affecting retes of Event 3: Competition Induced Death */    
    n = Sp * (Table->No_of_EVENTS-1) + k;   /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Vector;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Vector);
    #endif

    lambda_old   = Pa->rToI[n];  
    Pa->rToI[n] += Delta_Vector[i];
    Pa->rate[n] += f * Table->Competition_Induced_Death[Sp_0][j] / (double)Table->K_R;
    lambda_new   = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp, k, time_current, 
                                        lambda_old, lambda_new );                   
    #endif
  } 
  
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


