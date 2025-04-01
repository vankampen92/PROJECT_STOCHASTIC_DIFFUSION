#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, j, n, k;
  
  int Sp, R_0, R_1; 
  
  double K_R, y_S, m_0, D_0, D_1; 
  
  double Gam;           /* Conjugation/Encountering Rate in a Local Population  */
  double Dea;           /* Competition Induced Percapita Mortality Probability */
  double Eps;           /* Segregation Error at cell division */
  double Xhi;           /* Plasmid Transmission Probability */

  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites free from all species                   */
 
  Gam = Table->Lambda_R_1;  /* Conjugation/Encountering Rate in a Local Population  */
  Dea = Table->p_2;         /* Competition Induced Percapita Mortality Probability */
  Eps = Table->p_1;         /* Segregation Error at cell division */
  Xhi = Table->Chi_C_0;     /* Plasmid Transmission Probability */
 
  assert( Table->LOCAL_STATE_VARIABLES == Table->No_of_RESOURCES );

  for (j=0; j<Table->No_of_CELLS; j++) {
    
    y_S = Local_Population_Resources(j, y, Table);

    //assert(y_S <= K_R);

    m_0 = K_R-y_S;

    R_1   = j*Table->LOCAL_STATE_VARIABLES + 1;
    R_0   = j*Table->LOCAL_STATE_VARIABLES + 0;
     
    D_0   = Gam * Dea * y[R_1]/K_R;
    D_1   = Gam * Dea * y[R_0]/K_R; 

    /* Plasmid-free bacteria reproduce without segregation error:   Table->Beta_AP[0] * m_0/K_R *y[R]
       Plasmid-free bacteria reproduce with segregation error:      Table->Beta_AP[1] * m_0/K_R *y[R] * (1 - Eps)
    */
    dydt[R_0]  = Table->Beta_AP[0] * m_0/K_R *y[R_0] -Table->Delta_AP[0] *y[R_0]  -D_0 *y[R_0] -Gam * Xhi * y[R_0]/K_R * y[R_1] +Table->Beta_AP[1]*Eps * m_0/K_R *y[R_1];
      
    dydt[R_1]  = Table->Beta_AP[1]*(1.0 -Eps) * m_0/K_R *y[R_1] -Table->Delta_AP[1] *y[R_1] -D_1 *y[R_1] +Gam * Xhi * y[R_0]/K_R * y[R_1];
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


