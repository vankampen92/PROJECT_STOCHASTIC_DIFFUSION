#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

double Out_Mutation_Per_Capita(Parameter_Table * Table, int k, int j);
double In_Mutation(Parameter_Table * Table, int k, int j, const double * y); 

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j, k;
  int Sp; 
  double K_R, y_S, m_0; 

  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites free from all species                   */
 
  assert( Table->LOCAL_STATE_VARIABLES == 2 * Table->No_of_RESOURCES );

  for (j=0; j<Table->No_of_CELLS; j++) {
    
    y_S = Local_Population_Resources(j, y, Table);

    assert(y_S <= K_R);

    m_0 = K_R-y_S; 
    
    for(k = 0; k<Sp; k++) {
      
      RP  = j*Table->LOCAL_STATE_VARIABLES + 2*k;
      R   = j*Table->LOCAL_STATE_VARIABLES + 2*k+1;
    
      /* Lambda_R_0 represents an external passive
         source of propagules arriving in the local
         patches (at the same rate for all Sp, Lambda_R_0 ).
      */
      dydt[RP] = Table->Lambda_R_0 -Table->Delta_RP[k] *y[RP] + Table->Beta_AP[k] *y[R] -Table->Eta_RP[k] *m_0/K_R *y[RP];             
      dydt[R]  = -Table->Delta_AP[k] *y[R]  +Table->Eta_RP[k] *m_0/K_R *y[RP];
    }
  }

  if( Table->p_1 > 0.0) { /* Mutation along the trait axis */
    for (j=0; j<Table->No_of_CELLS; j++) {  
      for (k=0; k<Sp; k++) { 
        RP  = j*Table->LOCAL_STATE_VARIABLES + 2*k;
        R   = j*Table->LOCAL_STATE_VARIABLES + 2*k+1;

        dydt[RP] += 0.0;
        dydt[R]  += In_Mutation(Table, k, j, y) - Out_Mutation_Per_Capita(Table, k, j) * y[RP];
      }
    }
  }

  if( Table->No_of_CELLS > 1) { /* Physical Diffusion across Local Populations */  
    n= 0; 
    for (j=0; j<Table->No_of_CELLS; j++) { 
      
      for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++) { 
	      dydt[n] += In_Mu(Table, n, i, j, y) - Out_Mu_Per_Capita(Table, i, j) * y[n];
	      n++;
      }
    }
  }
  
  return GSL_SUCCESS;
}

double Out_Mutation_Per_Capita(Parameter_Table * Table, int k, int j)
{
  /* Input:
      . Table
      . k, Index of the Species 
      . j, Index of the Patch 
  
     Output:
      . Out_Mutation;
  */
  double Out_Mutation; 

  int i;
  int Sp; 
  double K_R, y_S, m_0;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
 
  Sp  = Table->No_of_RESOURCES;
  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites free from all species                   */

  double * y = Table->Vector_Model_Variables; 

  y_S = Local_Population_Resources(j, y, Table);
  m_0 = K_R-y_S;
  
  Out_Mutation = 2.0*Table->p_1 * m_0/K_R * Table->Eta_RP[k];

  return (Out_Mutation);
}

double In_Mutation(Parameter_Table * Table, int k, int j, const double * y)
{
  /* Input:
      . Table
      . k, Index of the Speciesz 
      . j, Index of the Patch 
  
     Output:
      . In_Mutation_Value;
  */
  double In_Mutation_Value; 

  int RP_Down, RP_Up; 
  int i;
  int Sp; 
  double K_R, y_S, m_0;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
 
  Sp  = Table->No_of_RESOURCES;
  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites free from all species                   */

  y_S = Local_Population_Resources(j, y, Table);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
  m_0 = K_R-y_S;

  i   = k-1; 
  RP_Down  = j*Table->LOCAL_STATE_VARIABLES + 2*i;
 
  i   = k+1; 
  RP_Up  = j*Table->LOCAL_STATE_VARIABLES + 2*i;

  if( k > 0 && k < Sp-1 )   
    In_Mutation_Value = 2.0*Table->p_1 * m_0/K_R  * (1.0/2.0* Table->Eta_RP[k-1]*y[RP_Down] + 1.0/2.0* Table->Eta_RP[k+1]*y[RP_Up]);
  else if( k == 0 )
    In_Mutation_Value = 2.0*Table->p_1 * m_0/K_R  * (Table->Eta_RP[k+1]*y[RP_Up]); 
  else if( k == Sp-1 )
    In_Mutation_Value = 2.0*Table->p_1 * m_0/K_R  * (Table->Eta_RP[k-1]*y[RP_Down]);
  else {
    /* Protection: just in case */
    printf("The index k of the strain/sp is out ouf range... ( 0 <= k <= Sp-1 ), but ( 0 <= k = %d <= Sp-1 )\n", 
            k);
    printf("The program will safely exit\n");
    exit(0);
  }

  return (In_Mutation_Value);
}

