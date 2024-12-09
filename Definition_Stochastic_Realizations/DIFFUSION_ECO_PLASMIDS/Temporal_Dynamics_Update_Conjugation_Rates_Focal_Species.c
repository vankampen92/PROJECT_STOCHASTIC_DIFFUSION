#include <MODEL.h>

#define TOLERANCE 1.0E-9

/* This functions allocate, initialize and free a number of local communities,
   which make up our total patch system or metapopulation */
extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Temporal_Dynamics_Update_Conjugation_Rates_Focal_Species ( Community * Pa, int x, Parameter_Table * Table, 
                                                                Stochastic_Rate * Rate, 
                                                                int Type_of_Event, int Type_Event, 
                                                                bool * bool_Next_Time, double * delta_rate )
{
   /* Changing rates across species/types as a result of a change in the amount of empty space (m_0). 
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
  double f, Delta_Vector; 
  int i, k, x, n, Sp, Sp_0;
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
  
  /* Conjugation Rates in Recipient List of the Focal Species */
  for( j = 0; j < Table->Recipient_List_Indeces[Sp_0][Table->No_of_RESOURCES]; j++) { 
    
    Sp = Table->Recipient_List_Indeces[Sp_0][j]; /* Strain ID of a recipient of species Sp_0 */
    assert( Sp != Sp_0); 

    Delta_Vector = f * P->Local_Strain_Population[Sp_0]->Gamma; / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
    
    n = Table->No_of_Event_Conjugation_Pair[Sp_0][Sp];  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Vector;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Vector);
    #endif

    lambda_old   = Pa->rToI[n];  
    Pa->rToI[n] += Delta_Vector[i];
    lambda_new   = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp_0, Sp, time_current, 
                                        lambda_old, lambda_new );                   
    #endif
  } 
  /* Conjugation Rates in Donor List of the Focal Species */
  for( j = 0; j < Table->Donor_List_Indeces[Sp_0][Table->No_of_RESOURCES]; j++) { 
    
    Sp = Table->Donor_List_Indeces[Sp_0][j]; /* Strain ID of a recipient of species Sp_0 */
    assert( Sp != Sp_0); 

    Delta_Vector = f * P->Local_Strain_Population[Sp]->Gamma; / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
           /* Delta Changes in conjugation rate of given conjugation pair  */

    n = Table->No_of_Event_Conjugation_Pair[Sp][Sp_0];  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Vector;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Vector);
    #endif

    lambda_old     = Pa->rToI[n];  
      Pa->rToI[n] += Delta_Vector[i];
    lambda_new     = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp_0, Sp, time_current, 
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


