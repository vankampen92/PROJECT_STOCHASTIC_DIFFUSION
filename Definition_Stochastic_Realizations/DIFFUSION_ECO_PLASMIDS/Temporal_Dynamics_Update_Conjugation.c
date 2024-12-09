#include <MODEL.h>

#define TOLERANCE 1.0E-9

extern gsl_rng * r; /* Global generator defined in main.c */
#define RANDOM gsl_rng_uniform_pos(r)

void Temporal_Dynamics_Update_Conjugation ( Community * Pa, int x, Parameter_Table * Table, int * Patch, 
                                            Stochastic_Rate * Rate, int Type_of_Event,  
                                            bool * bool_Next_Time, double * delta_rate )
{
   /* Changing rates across recipient/transconjugant and its donors and recipients as a result 
      of a successful conjugation event. 
    */
    /* Input: 
      . Pa, focal patch 
      . x, index of focal patch
      . Table
      .  Patch is an array containing the two patches involved in a movement event. 
                    Patch[0] is the patch sending the individual
		                Patch[1] is the patch receiving the individual.
                    Patch[2] is the strain ID of the RECIPIENT species
                    Patch[3] is the strain ID of the TRANSCONJUGANT species
                    Patch[4] is 0 if no update is required and 1 otherwise. 
      . Rate
      . Type_of_Event

      Output: 
      . Delta_Rate.
      . bool_Next_Time. 
    */ 
  double nR, nT; 
  double B0;
  double B1;
  int Adjacent_List[4] = {0, 2, 4, 5};
  double Delta_Vector[4];
  double Delta_Value; 

  int i, k, x, n, N, Sp, Sp_Recip, Sp_Trans;

  int Index_Leaf_x, Index_Node_x;
  int Index_Leaf_y, Index_Node_y; 
  double time_old, lambda_old, lambda_new, Next_Time;
  treenode * Leaf;
  treenode * Node;   

  double time_current = Rate->Stochastic_Time;
  double Delta_Rate = * delta_rate; 
  
  N = Type_of_Event - Table->No_of_RESOURCES * (Table->No_of_EVENTS-1); 
  assert( N < Table->No_of_CONJUGATION_EVENTS );
  
  /* Sp: strain ID of the RECIPIENT species */
  Sp_Recip = Patch[2];   
  Sp       = Patch[2];  
  assert( Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[N][1] == Sp );

  B0  = Table->Beta_AP[Sp] * (1.0 - Table->p_1);
  B1  = Table->Beta_AP[Sp] * Table->p_1;
  Delta_Vector[0] = -Table->Mu_RP[Sp];
  Delta_Vector[1] = -Table->Delta_AP[Sp];
  Delta_Vector[2] = -B0 / (double)Table->K_R * (double)Pa->m_0;
  Delta_Vector[3] = -B1 / (double)Table->K_R * (double)Pa->m_0;

  for(i=0; i<4; i++) {
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

  /* Sp: the strain ID of the TRANSCONJUGANT species */
  Sp_Trans = Patch[3]; 
  Sp       = Patch[3];
  B0  = Table->Beta_AP[Sp] * (1.0 - Table->p_1);
  B1  = Table->Beta_AP[Sp] * Table->p_1;
  Delta_Vector[0] = +Table->Mu_RP[Sp];
  Delta_Vector[1] = +Table->Delta_AP[Sp];
  Delta_Vector[2] = +B0 / (double)Table->K_R * (double)Pa->m_0;
  Delta_Vector[3] = +B1 / (double)Table->K_R * (double)Pa->m_0;
  
  for(i=0; i<4; i++) {
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

  /* Changes in the conjugation rates of the both Recipient and Transconjugant with, in turn, all its 
     respective recipient as well as in the donors of them */

  /* Changes in Conjugation Rates in Recipient List of the Focal RECIPIENT Species */
  for( j = 0; j < Table->Recipient_List_Indeces[Sp_Recip][Table->No_of_RESOURCES]; j++) { 
  
    Sp = Table->Recipient_List_Indeces[Sp_Recip][j]; /* Strain ID of a recipient of species Sp_Recip */
    assert( Sp != Sp_Recip); 

    if(Sp != Sp_Trans) 
      Delta_Value = - P->Local_Strain_Population[Sp_Recip]->Gamma; / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
    else {
      nT = (double)Pa->Local_Strain_Population[Sp_Trans]->n;
      nR = (double)Pa->Local_Strain_Population[Sp_Recip]->n;  
      Delta_Value = P->Local_Strain_Population[Sp_Recip]->Gamma; / (double)Table->K_R * ( nR - nT + 1.0 );
    } 

    n = Table->No_of_Event_Conjugation_Pair[Sp_Recip][Sp];  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Value;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Value);
    #endif

    lambda_old   = Pa->rToI[n];  
    Pa->rToI[n] += Delta_Value;
    lambda_new   = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp_Recip, Sp, time_current, 
                                        lambda_old, lambda_new );                   
    #endif
  }

  /* Changes in Conjugation Rates in Donor List of the Focal RECIPIENT Species */
  for( j = 0; j < Table->Donor_List_Indeces[Sp_Recip][Table->No_of_RESOURCES]; j++) { 
    
    Sp = Table->Donor_List_Indeces[Sp_Recip][j];     /* Strain ID of a Donor List of species Sp_Recip */
    assert( Sp != Sp_Recip); 

    f(Sp != Sp_Trans) 
      Delta_Value = - P->Local_Strain_Population[Sp]->Gamma; / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
    else {
      nT = (double)Pa->Local_Strain_Population[Sp_Trans]->n;
      nR = (double)Pa->Local_Strain_Population[Sp_Recip]->n;  
      Delta_Value = P->Local_Strain_Population[Sp_Trans]->Gamma; / (double)Table->K_R * ( nR - nT + 1.0 );
    }

    n = Table->No_of_Event_Conjugation_Pair[Sp][Sp_Recip];  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Value;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Value);
    #endif

    lambda_old     = Pa->rToI[n];  
    Pa->rToI[n] += Delta_Value;
    lambda_new     = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp, Sp_Recip, time_current, 
                                        lambda_old, lambda_new );                   
    #endif
  }

  /* Changes in Conjugation Rates in Recipient List of the Focal TRANCONJUGANT Species */
  for( j = 0; j < Table->Recipient_List_Indeces[Sp_Trans][Table->No_of_RESOURCES]; j++) { 
    
    Sp = Table->Recipient_List_Indeces[Sp_Trans][j]; /* Strain ID of a recipient of species Sp_Trans */
    assert( Sp != Sp_Trans); 

    if( Sp != Sp_Recip)
      Delta_Value = P->Local_Strain_Population[Sp_Trans]->Gamma; / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
    else {
      nT = (double)Pa->Local_Strain_Population[Sp_Trans]->n;
      nR = (double)Pa->Local_Strain_Population[Sp_Recip]->n;  
      Delta_Value = P->Local_Strain_Population[Sp_Trans]->Gamma; / (double)Table->K_R * ( nR - nT + 1.0 );      
    } 

    n = Table->No_of_Event_Conjugation_Pair[Sp_Trans][Sp];  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Value;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Value);
    #endif

    lambda_old   = Pa->rToI[n];  
    Pa->rToI[n] += Delta_Value;
    lambda_new   = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp_Trans, Sp, time_current, 
                                        lambda_old, lambda_new );                   
    #endif
  } 

  /* Changes in Conjugation Rates in Donor List of the Focal TRANSCONJUGANT Species */
  for( j = 0; j < Table->Donor_List_Indeces[Sp_Trans][Table->No_of_RESOURCES]; j++) { 
    
    Sp = Table->Donor_List_Indeces[Sp_Trans][j]; /* Strain ID of a donor list of species Sp_Trans */
    assert( Sp != Sp_Trans); 

     if( Sp != Sp_Recip)
      Delta_Value = P->Local_Strain_Population[Sp]->Gamma; / (double)Table->K_R * (double)Pa->Local_Strain_Population[Sp]->n;
     else {
      nT = (double)Pa->Local_Strain_Population[Sp_Trans]->n;
      nR = (double)Pa->Local_Strain_Population[Sp_Recip]->n;  
      Delta_Value = P->Local_Strain_Population[Sp_Recip]->Gamma; / (double)Table->K_R * ( nR - nT + 1.0 );
     }

    n = Table->No_of_Event_Conjugation_Pair[Sp][Sp_Trans];  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */

    Delta_Rate += Delta_Value;

    #if defined BINARY_TREE_SUPER_OPTIMIZATION
      Index_Leaf_x = Table->TOTAL_No_of_EVENTS * x + n;
      Leaf = Table->Leaves[Index_Leaf_x];
      sum_Delta_upto_Root(Table->Treeroot, Leaf, Delta_Value);
    #endif

    lambda_old     = Pa->rToI[n];  
    Pa->rToI[n] += Delta_Value;
    lambda_new     = Pa->rToI[n];
    
    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      Priority_Queu_Super_Optimization( Pa, x, Table, n, Sp, Sp_Trans, time_current, 
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

  free(Delta_Vector);    
}


