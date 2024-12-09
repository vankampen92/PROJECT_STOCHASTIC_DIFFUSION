#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j, k;
  int Sp; 
  double K_R, y_S, m_0; 
  double DE ;           /* Competition Induced Percapita Mortality Rate */
  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; /* Total Number of Local Sites in each Local Population                      */
  y_S = 0.0;                /* Total Population (Local Sites Occupied by all species from k=0 to k=Sp-1) */
  m_0 = 0.0;                /* Total Number of Local Empty Sites free from all species                   */
 
  assert( Table->LOCAL_STATE_VARIABLES == Table->No_of_RESOURCES );

  for (j=0; j<Table->No_of_CELLS; j++) {
    
    y_S = Local_Population_Resources(j, y, Table);

    assert(y_S <= K_R);

    m_0 = K_R-y_S;

    
    /* Segregation Total Error (across strains IDs of the same bacterial type) */
    double * ER = (double *)calloc(Table->No_of_STRAINS, sizeof(double));    
    int * Ix    = (int *)calloc(Table->No_of_RESOURCES, sizeof(int));     
    /* Total Rate of creation of new transconjugants */
    double * CT_Gain = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));    
    /* Total Rate of loss of recipients */
    double * CR_Loss = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));

    Transconjugation_Gain_and_Loss_Total_Rates( y, CT_Gain, CR_Loss, Table );

    k = 0; 
    for(i = 0; i < Table->No_of_STRAINS; i++) {
      
      ER[i] = 0.0; 
      for( n = 0; n < Table->n[i]; n++ ) { 
        R   = j*Table->LOCAL_STATE_VARIABLES + k;
 
        ER[i] += Table->Beta_AP[k] * (Table->p*1) * m_0/K_R *y[R];
        
        Ix[k] = i;          
        k++;
      }
    } 
    
    assert (k == Sp);

    DE = 0.0;   
    for(k = 0; k<Sp; k++) {

      R   = j*Table->LOCAL_STATE_VARIABLES + k;
      DE  = Competition_Induced_Death_Based_Calculation(Table, k, j);

      /* Plasmid-free bacteria reproduce without segregation error: 
           
           For the strains IDs corresponding to plasmid-free profiles, the equation would look like:
           
           dydt[R] = Table->Beta_AP[k] * (1.0 - Table->p*1)*m_0/K_R *y[R] + Table->Beta_AP[k]*(Table->p*1) + m_0/K_R *y[R] + ... 
           
           where the last term stems from a term in the ER[Ix[k]] accounting from replenishment of strain-free IDs
           through segregation error, which, in case of plasmid free Strain IDs, amounts to a total increasing rate 
           (by summing the two terms) through cell division without segregation error:
           
           Table->Beta_AP[k] * m_0/K_R *y[R]

           (of course, in plasmid-free bacteria there are not plasmids to segregate during cell division) 
      */

      if ( k == Table->n_0[Ix[k]] )
        dydt[R]  = Table->Beta_AP[k] * (1.0 -Table->p*1) * m_0/K_R *y[R] +ER[Ix[k]] -Table->Delta_AP[k] *y[R] -DE *y[R] -CR_Loss[R];
      else
        /* Plasmid-nonfree bacteria reproduce with segregation error */
        dydt[R]  = Table->Beta_AP[k] * (1.0 -Table->p*1) * m_0/K_R *y[R] -Table->Delta_AP[k] *y[R] -DE *y[R]  +CT_Gain[R] -CR_Loss[R];       
    }

    free(Ix); 
    free(ER);
    free(CR_Loss);
    free(CT_Gain);
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


