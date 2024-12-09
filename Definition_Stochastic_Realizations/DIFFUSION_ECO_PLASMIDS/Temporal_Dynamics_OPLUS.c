/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2010 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

double Competition_Induced_Mortality_Calculation(Parameter_Table * , int , int );
double  Conjugation_Percapita_Rate_Calculation(Parameter_Table * , int , int );

void Temporal_Dynamics(Community ** My_Community, Parameter_Table * Table, Stochastic_Rate * Rate)
{
  /* This function calculates the rates of all possible events from scratch, across and within cells 
     given certain system confituration as defined in My_Community 
  */
  /* For this function, the calling hierarchy up to this level is: 
     main() (main.c) 
        --->  M_O_D_E_L___S_T_O( &Table ); (MODEL_STO.c) 
            --->  S_T_O_C_H_A_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S (...); (Stochastic_Time_Dynamics.c)
                ---> Temporal_Dynamics(...) (./Definition_Stochastic_Realizations/DIFFUSION_ECO_PLASMIDS/Temporal_Dynamics.c) 
  */
  int i,j,k,n, N, Sp;
  Community * P;
  int No_of_CELLS;
  int GRAND_No_of_EVENTS;
  double OutMigration;
  double K_R, m_0, y_S; 
  double Competition_Induced_Percapita_Mortality_Rate; 
  double Conjugation_Percapita_Rate; 

  Parameter_Model * pa  = Table->P;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites (free from all species)                 */
  
  P = My_Community[0];  /* P could be used as a pointer to the zero-th to be incremented 
			                     if necessary (not used like that in this implementation) 
		                    */
  No_of_CELLS             = pa->No_of_CELLS;
  Sp                      = pa->No_of_RESOURCES; 
    
  double * Y              = Table->Vector_Model_Variables; 
  
  Immigration_Preassure_on_Focal_Patch_Initialization( My_Community, pa );

  Rate->max_Probability = 0.0;
  Rate->Total_Rate      = 0.0;
  
  GRAND_No_of_EVENTS = 0;
  for(i=0; i<No_of_CELLS; i++){

    P = My_Community[i];
    y_S = Local_Population_Resources(i, Y, Table);
    m_0 = K_R-y_S;                                /* Total empty microsites (local sites) in the i-th cell */

    assert( (double)P->m_0 == m_0 );

    P->ratePatch = 0; 
    n = 0;
    for(k=0; k<Sp; k++) {

      assert(Table->R == 0);
      R  = k + Table->R; 
      /* Strain-profile type (from 0 to Sp-1). Notice that Table->R should be zero!!! */
      assert(R == k);

      /* 0: Bacteria Out-Migration (R --> R-1) and some other patch gains one */ 
      OutMigration = P->Total_Per_Capita_Out_Migration_Rate[R];                    
      assert( 4 * Table->Mu ==  OutMigration );                 /* Movements on a squared lattice */
          
      P->rate[n] = OutMigration;       
      P->rToI[n] = OutMigration * (double)P->n[R]; 
      P->ratePatch += P->rToI[n];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif 
      n++;

      /* 1: Bacteria External Immigration event */
      P->rate[n] = Table->Lambda_R_0;  
      P->rToI[n] = Table->Lambda_R_0 * m_0; 
      P->ratePatch += P->rToI[n]; 

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;

      /* 2: Bacteria Death  (R --> R-1) */
      P->rate[n] = Table->Delta_AP[k];   
      P->rToI[n] = Table->Delta_AP[k] * (double)P->n[R];
      P->ratePatch += P->rToI[n];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;

      /* 3: Competition induced mortality on Species R  (R --> R-1) */
      Competition_Induced_Percapita_Mortality_Rate = Competition_Induced_Death_Based_Calculation(Table, R, i);
      P->rate[n] = Competition_Induced_Percapita_Mortality_Rate;  
      P->rToI[n] = Competition_Induced_Percapita_Mortality_Rate * (double)P->n[R];
      P->ratePatch += P->rToI[n];
      
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;

      /* 4: Cell dividion (without error) */
      P->rate[n] = Table->Beta_AP[k] * (1.0 - Table->p_1) * m_0/K_R;     
      P->rToI[n] = Table->Beta_AP[k] * (1.0 - Table->p_1) * m_0/K_R * (double)P->n[R];
      P->ratePatch += P->rToI[n];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
      
      /* 5: Cell division (with error) */
      P->rate[n] = Table->Beta_AP[k] * (Table->p_1) * m_0/K_R;     
      P->rToI[n] = Table->Beta_AP[k] * (Table->p_1) * m_0/K_R * (double)P->n[R];
      P->ratePatch += P->rToI[n];

      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
      #endif
      #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
        Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
      #endif
      n++;
    }

    N = 0; 
    for(k=0; k<Sp; k++) {
      /* 6: Cell to Cell Conjugation Events */  
      /* Donor (k, p) + Recipient (j, q) --->  Donor (k, p) + Transconjugant (l, r) */
      R = k; 

      for(j=0; j<P->Local_Strain_Population[k]->Recipient_List[Sp]; j++) {
        
        /* Event: 6.j: Conjugation of the donor Strain ID = k (or R) with j-th recipient 
           (Strain ID = Strain_ID_1) in the recipent list of Strain ID = k (or R)
        */
        Strain_ID_1  = P->Local_Strain_Population[k]->Recipient_List[j];   

        /* Per capita rate: */
        P->rate[n] = P->Local_Strain_Population[k]->Gamma * (double)P->Local_Strain_Population[Strain_ID_1]->n / K_R;  
        /* Propensity (or total transition rate for the event 6.j): */  
        P->rToI[n]  = P->rate[n] * (double)P->n[R];
        P->ratePatch += P->rToI[n];

        Table->No_of_Event_Conjugation_Pair[k][Strain_ID_1] = n;  /* n >= 0  and n < Table->TOTAL_No_of_EVENTS */
        
        Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[N][0] = k;            /* DONOR     */
        Table->Event_Conjugation_Donor_Recipient_Pair_Strain_IDs[N][1] = Strain_ID_1;  /* RECIPIENT */

        #if defined BINARY_TREE_SUPER_OPTIMIZATION
          Table->Leaves[GRAND_No_of_EVENTS++]->value = P->rToI[n];
        #endif
        #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
          Table->T->Vector_of_Rates[GRAND_No_of_EVENTS++] = P->rToI[n];  
        #endif
        n++;

        N++; 
      }        
    }
    assert( N == Table->No_of_CONJUGATIONS_EVENTS );

    #if defined BINARY_TREE_OPTIMIZATION
      Table->Leaves[i]->value = P->ratePatch;
    #endif

    assert( n == Table->TOTAL_No_of_EVENTS );
    
    Rate->Total_Rate += P->ratePatch;
    Rate->max_Probability = MAX( Rate->max_Probability, P->ratePatch );
  }

  #if defined BINARY_TREE_SUPER_OPTIMIZATION
  assert(GRAND_No_of_EVENTS == Table->TOTAL_GRAND_No_of_EVENTS);
  #endif 
  
  #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
  assert(GRAND_No_of_EVENTS == Table->TOTAL_GRAND_No_of_EVENTS);
  #endif 

  if(Rate->Total_Rate <= 0.0){
      printf("\n");
      printf(" R is the total temporal rate of system configuration change\n");
      printf(" R = %g\n", Rate->Total_Rate );
      printf(" As R is zero or negative, no change is possible\n");
      printf(" R shouldn't be negative. If it is, there are definitely some errors in the code\n");
      printf("\n");
      if( Rate->Total_Rate < 0.0 ) exit(0);
  }
}

double Competition_Induced_Death_Based_Calculation(Parameter_Table * Table, 
                                                   int Strain_ID, int Patch)
{
  double Rate; 
  int i,j, Sp;
  int i_Str, j_Str, i_Pro, j_Pro;
  int Strain_ID_1; 
  double K_R; 

  Parameter_Model * pa  = Table->P;
  Community * P         = Table->Patch_System[Patch];

  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  Sp = pa->No_of_RESOURCES; 
  
  Rate = 0.0; 
  for(j=0; j<P->Local_Strain_Population[Strain_ID]->Competition_List[Sp]; j++) {
    Strain_ID_1  = P->Local_Strain_Population[Strain_ID]->Competition_List[j];

    Rate += Table->Competition_Induced_Death[Strain_ID][j]*(double)P->Local_Strain_Population[Strain_ID_1]->n / K_R;  
  }  
    
  return(Rate);
}
