#include <MODEL.h>

#define TOLERANCE 1.0E-9

/* This functions allocate, initialize and free a number of local communities,
   which make up our total patch system or metapopulation */
extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Temporal_Dynamics_Update_Empty_Space_Induced ( Community * Pa, int x, Parameter_Table * Table, 
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
  double f; 
  double B0;
  double B1;
  int Adjacent_List[3] = {1, 4, 5};

  int i, k, x, n, Sp, Sp_0;

  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double time_old, lambda_old, lambda_new, Next_Time;
  treenode * Leaf;
  treenode * Node;   

  double time_current = Rate->Stochastic_Time;
  double Delta_Rate = * delta_rate; 
  
  double * Delta_Vector = (double *)calloc(Table->No_of_EVENTS); 
     
  if(Type_Event == 0 || Type_Event == 2 || Type_Event == 3)
        f = ( 1.0);
  else if (Type_Event == 1 || Type_Event == 4)
        f = (-1.0);
  else {
        printf("Error in Temporal_Dynamics_Empty_Space_Induced.c\n");
        exit(0); 
  }

  Delta_Vector[0] = f *Table->Lambda_R_0;

  Sp_0  = Type_of_Event/(Table->No_of_EVENTS-1);

  for( Sp = 0; Sp < Table->No_of_RESOURCES; Sp++) 
    if( Sp != Sp_0) {

      B0  = Table->Beta_AP[Sp] * (1.0 - Table->p_1);
      B1  = Table->Beta_AP[Sp] * Table->p_1;

      Delta_Vector[1] = f * B0 / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
      Delta_Vector[2] = f * B1 / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
       
      for(i=0; i<3; i++) {
        k = Adjacent_List[i];

        Delta_Rate += Delta_Vector[i];

        #if defined BINARY_TREE_SUPER_OPTIMIZATION
          Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + Sp*(Table->No_of_EVENTS-1) + k;
          Leaf = Table->Leaves[Index_Leaf_x];
          sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Vector[i]);
        #endif

        lambda_old                                = Pa->rToI[Sp*(Table->No_of_EVENTS-1) + k];  
        Pa->rToI[Sp*(Table->No_of_EVENTS-1) + k] += Delta_Vector[i];
        lambda_new                                = Pa->rToI[Sp*(Table->No_of_EVENTS-1) + k];
    
        #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
          n = Sp*(Table->No_of_EVENTS-1) + k;  
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

  free(Delta_Vector);    
}


