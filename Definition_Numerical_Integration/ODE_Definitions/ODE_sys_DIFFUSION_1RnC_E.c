#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, m, j;
  int Sp; 
  double K_R;
  double A_T, A_F;  /* Total and Total Adult (fertile) Population */ 
  double RA_T;      /* Total Compound Population */
  double ARA_CC_0, ARA_CD_0, ARA_DC_0, ARA_DD_0;
  double ARA_CC_1, ARA_CD_1, ARA_DC_1, ARA_DD_1;
  double p_1, p_2, y_A; 
  
  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES;
  assert( Sp == 1 ); 

  int *    A = (int *)calloc(Table->N_E, sizeof(int));
  int *   RA = (int *)calloc(Table->N_E, sizeof(int));
  int ** ARA = (int **)calloc(Table->N_E, sizeof(int *));
  for(i=0; i<Table->N_E; i++) 
    ARA[i] = (int *)calloc(Table->N_E, sizeof(int));
  
  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  K_R = (double)Table->K_R; 
  p_1 = Table->p_1;
  p_2 = Table->p_2; 
  
  for (j=0; j<Table->No_of_CELLS; j++) {

    R      = j*Table->LOCAL_STATE_VARIABLES + Table->R;

    for(n=0; n<Table->N_E; n++) {
      A[n] = j*Table->LOCAL_STATE_VARIABLES + Table->A_P[n];
      RA[n]= j*Table->LOCAL_STATE_VARIABLES + Table->RA_P[n];
    }
    
    for(n=0; n<Table->N_E; n++)
      for(m=0; m<Table->N_E; m++)
	ARA[n][m] = j*Table->LOCAL_STATE_VARIABLES + Table->ARA_P[n][m];

    A_T = 0.0; A_F = 0.0; RA_T = 0.0; 
    for(i=0; i<Table->N_E; ++) {
      A_T  += y[A[i]];
      RA_T += y[RA[i]]; 
      if( i >= Table->i_0 ) A_F += y[A[i]];  /* Only Mature Population contributing to reproduction */
    }
    
    dydt[R] = Table->Lambda_R_0 *(K_R-y[R]) +Table->Beta_R *(K_R-y[R])/K_R *y[R] -Table->Delta_R_0 *y[R] -Table->Alpha_C_0 *y[R]/K_R * A_T;

    for(i=0; i<Table->N_E; ++) {
      
                      ARA_CD_0 = 0.0; ARA_DC_0 = 0.0; ARA_DD_0 = 0.0;
      ARA_CC_1 = 0.0; ARA_CD_1 = 0.0; ARA_DC_1 = 0.0; ARA_DD_1 = 0.0;
      for(n=0; n<Table->N_E; n++) {
	ARA_CD_0 += Table->Eta_C_0* p_1*(1.0-p_2)*       y[ARA[i][n]] ;
	ARA_DC_0 += Table->Eta_C_0* (1.0-p_1)*p_2*       y[ARA[n][i]] ;
	ARA_DD_0 += Table->Eta_C_0* (1.0-p_1)*(1.0-p_2)* y[ARA[i][n]] ;

	if( i >= Table->k_E) {  
	  ARA_CC_1 += Table->Eta_C_0* p_1*p_2* (y[ARA[i-Table->k_E][n]] + y[ARA[n][i-Table->k_E]]);
	  ARA_DD_1 += Table->Eta_C_0* (1.0-p_1)*(1.0-p_2)* y[ARA[n][i-Table->k_E]];
	}
	if( i >= 2*Table->k_E) {
	  ARA_DC_1 += Table->Eta_C_0* (1.0-p_1)*p_2* y[ARA[i-2*Table->k_E][n]];
	  ARA_CD_1 += Table->Eta_C_0* p_1*(1.0-p_2)* y[ARA[n][i-2*Table->k_E]];
	}
      }

      if(i == 0) 
	dydt[A[0]] = Table->Beta_C*Table->f*A_F +Table->Lambda_C_0 +Table->Theta_C*y[A[1]] -Table->Delta_C_0 *y[A[0]] -Table->Alpha_C_0 *y[R]/K_R *y[A[0]] -Table->Chi_C_0 * RA_T/K_R *y[A[0]] + ARA_CD_0 + ARA_DC_0 + ARA_DD_0;
      else if (i < Table->k_E) 
	dydt[A[i]] = Table->Lambda_C_0 + Table->Theta_C*y[A[i+1]] -Table->Theta_C*y[A[i]] -Table->Delta_C_0 *y[A[i]] -Table->Alpha_C_0 *y[R]/K_R *y[A[i]] -Table->Chi_C_0 * RA_T/K_R *y[A[i]] + ARA_CD_0 + ARA_DC_0 + ARA_DD_0;
      else if (i < 2*Table->k_E)
	dydt[A[i]] = Table->Lambda_C_0 + Table->Theta_C*y[A[i+1]] -Table->Theta_C*y[A[i]] -Table->Delta_C_0 *y[A[i]] -Table->Alpha_C_0 *y[R]/K_R *y[A[i]] -Table->Chi_C_0 * RA_T/K_R *y[A[i]] + ARA_CD_0 + ARA_DC_0 + ARA_DD_0 + ARA_CC_1 + ARA_DD_1;
      else if (i < (Table->N_E-2))
	dydt[A[i]] = Table->Lambda_C_0 + Table->Theta_C*y[A[i+1]] -Table->Theta_C*y[A[i]] -Table->Delta_C_0 *y[A[i]] -Table->Alpha_C_0 *y[R]/K_R *y[A[i]] -Table->Chi_C_0 * RA_T/K_R *y[A[i]] + Table->Nu_C_0 * y[A[i-2*Table->k_E]] + ARA_CD_0 + ARA_DC_0 + ARA_DD_0 + ARA_CC_1 + ARA_CD_1 + ARA_DC_1 + ARA_DD_1;
      else {
	ARA_CC_1 = 0.0; ARA_CD_1 = 0.0; ARA_DC_1 = 0.0; ARA_DD_1 = 0.0;
	for(n=0; n<Table->N_E; n++) { 
	  for(m=0; m<=Table->k_E; m++) { 
	    ARA_CC_1 += Table->Eta_C_0* p_1*p_2* (y[ARA[i-m][n]] + y[ARA[n][i-m]]);
	    ARA_DD_1 += Table->Eta_C_0* (1.0-p_1)*(1.0-p_2)* y[ARA[n][i-m]];
	  }
	  for(m=0; m<=2*Table->k_E; m++) { 
	    ARA_DC_1 += Table->Eta_C_0* (1.0-p_1)*p_2* y[ARA[i-m][n]];
	    ARA_CD_1 += Table->Eta_C_0* p_1*(1.0-p_2)* y[ARA[n][i-m]];
	  }
	}
      
	y_RA = 0.0; 
	for(m=0; m<=2*Table->k_E; m++) y_RA += y[RA[i-m]];  
	  
	dydt[A[i]] = Table->Lambda_C_0 -Table->Theta_C*y[A[i]] -Table->Delta_C_0 *y[A[i]] -Table->Alpha_C_0 *y[R]/K_R *y[A[i]] -Table->Chi_C_0 * RA_T/K_R *y[A[i]] + Table->Nu_C_0 *y_RA + ARA_CD_0 + ARA_DC_0 + ARA_DD_0 + ARA_CC_1 + ARA_CD_1 + ARA_DC_1 + ARA_DD_1;
		
      }   
    }

    for(i=0; i<Table->N_E; ++)
      dydt[RA[i]] = Table->Alpha_C_0 *y[R]/K_R *y[A[i]] -Table->Nu_C_0 *y[RA[i]] -Table->Chi_C_0 * y[RA[i]]/K_R *A_T; 
    
    for(n=0; n<Table->N_E; n++)
      for(m=0; m<Table->N_E; m++) 
	dydt[ARA[n][m]] = Table->Chi_C_0 * y[RA[m]]/K_R *y[A[n]] - Table->Eta_C_0 * y[ARA[n][m]]; 
    
  }
  
  n= 0; 
  for (j=0; j<Table->No_of_CELLS; j++) { 
    
    for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++) { 
      dydt[n] += In_Mu(Table, n, i, j, y) - Out_Mu_Per_Capita(Table, i, j) * y[n];;
      n++;
    }
  }
    
  for(i=0; i<Table->N_E; i++) free(ARA[i]);
  free(ARA); free(A); free(RA); 
  
  return GSL_SUCCESS;
}
