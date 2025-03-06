double Gam = Table->Lambda_R_1;  /* Conjugation/Encountering Rate in a Local Population  */
double Dea = Table->p_2;         /* Competition Induced Percapita Mortality Probability */
double Eps = Table->p_1;         /* Segregation Error at cell division */
double Xhi = Table->Chi_C_0;     /* Plasmid Transmission Probability */

assert( Table->LOCAL_STATE_VARIABLES == Table->No_of_RESOURCES );
assert( Table->No_of_CELLS == 1 );

double x;
double K_R = (double)Table->K_R;

/* F_0( y; parameters) */
  x = Table->Beta_AP[0] - Table->Delta_AP[0]  
    - 2.0 * Table->Beta_AP[0] / K_R * y[0] 
    - Table->Beta_AP[0] * y[1]/K_R
    - Gam * Dea * 0.5 * y[1]/K_R 
    - Gam * Xhi * y[1]/K_R 
    - Table->Beta_AP[1] * Eps * y[1]/K_R;
  gsl_matrix_set(m, 0, 0, x);

  x = Table->Beta_AP[1]*Eps - 2.0 * Table->Beta_AP[1]*Eps * y[1]/K_R;
    - Table->Beta_AP[0] * y[0]/K_R 
    - 0.5 * Gam * Dea * y[0]/K_R 
    - Gam * Xhi * y[0]/K_R
    - Table->Beta_AP[1]*Eps * y[0]/K_R; 
  gsl_matrix_set(m, 0, 1, x);

  /* F_1( y; parameters) */
  x = -Table->Beta_AP[1]*(1.0 -Eps) * y[1]/K_R 
    - Gam * Dea * 0.5 *y[1]/K_R 
    + Gam * Xhi * y[1]/K_R;
  gsl_matrix_set(m, 1, 0, x);

  x = Table->Beta_AP[1]*(1.0 -Eps) - Table->Delta_AP[1]
    - 2.0 * Table->Beta_AP[1]*(1.0 -Eps) * y[1]/K_R
    - Table->Beta_AP[1]*(1.0 -Eps) * y[0]/K_R 
    - 0.5 * Gam * Dea * y[0]/K_R  
    + Gam * Xhi * y[0] /K_R;
  gsl_matrix_set(m, 1, 1, x);

