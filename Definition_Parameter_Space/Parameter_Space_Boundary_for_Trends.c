#include <MODEL.h>

void Parameter_Space_Boundary_for_Trends(int j, double *value_0, double *value_1)
{
  switch(j)
  {
    case 0: 
    (*value_0) = 0.0;        (*value_1) = 600.0;
    break;
    case 1: 
    (*value_0) = 0.0;        (*value_1) = 1.0;
    break;
    case 2: 
    (*value_0) = 0.0;        (*value_1) = 50000.0;
    break;
    case 3: 
    (*value_0) = (1.0/70.0); (*value_1) = (1.0/20.0);
    break;
    case 4: 
    (*value_0) = (1.0/70.0); (*value_1) = (1.0/20.0);
    break;
    case 5: 
    (*value_0) = (6.0);      (*value_1) = (24.0);
    break;
    case 6: 
    (*value_0) = (1.0/20.0); (*value_1) = (1.0/5.0);
    break;
    case 7: 
    (*value_0) =  0.0;       (*value_1) =  99.0;
    break;
    case 8: 
    (*value_0) =  0.0;       (*value_1) =  600.0;
    break;
    /* case 9: /\* p_YX   *\//\* Female-to-Male transmission probability T(X--->Y)  *\//\*   9 *\/ */
    /* (*value_0) = 0.0;        (*value_1) = 1.0; */
    /* break; */
    /* case 10: /\* Sigma_0   *\//\* Transitional rate into the sexual worker stage *\//\*  10 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (1.0/2.0); */
    /* break; */
    /* case 11: /\* Sigma_0_r *\//\* Reversal rate back into the regural female stage*\//\* 11 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (1.0/2.0); */
    /* break; */
    /* case 12: /\* Eta     *\//\* Contact rate relative increase for sexual workers *\//\* 12 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (9.0); */
    /* break; */
    /* case 13: /\* F_X     *\//\* Total no of females entering sexual age per year  *\//\* 13 *\/ */
    /* (*value_0) =   0.0;      (*value_1) = 50000.0; */
    /* break; */
    /* case 14: /\* n_Mu   *\//\* No of latent subclasses during quasi-latent phase  *\//\* 14 *\/ */
    /* (*value_0) = ( 1.0);     (*value_1) = (10.0); */
    /* break; */
    /* case 15: /\* Xhi   *\//\* Transmision probability increase during acute phase *\//\* 15 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (9.0); */
    /* break; */
    /* case 16: /\* f_W   *\//\* Fraction of male sexual contacts with females sexual workers *\//\* 16 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (0.5); */
    /* break; */
    /* case 17: /\* f_0   *\//\* Fraction of male sexual contacts with young females *\/ /\* 17 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (0.5); */
    /* break; */
    /* case 18: /\* Sigma_1   *\//\* Transitional rate into the sexual worker stage *\/  /\* 18 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (1.0/2.0); */
    /* break; */
    /* case 19: /\* Sigma_1_r *\//\* Reversal rate back into the regural female stage *\//\* 19 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (1.0/2.0); */
    /* break; */
    /* case 20: /\* Alpha     *\//\* Aging: Transition rate into adult stage *\/         /\* 20 *\/ */
    /* (*value_0) = (1.0/15.0); (*value_1) = (1.0/10.0); */
    /* break; */
    /* case 21: /\* r_F     *\//\* Factor determining fertile population fraction  *\/   /\* 21 *\/ */
    /* (*value_0) = (0.000001); (*value_1) = (1.0); */
    /* break; */
    /* case 22: /\* D_Sigma   *\//\* Parameter introducing time-dependence in Sigma_0 and Sigma_1 *\//\* 22 *\/ */
    /* (*value_0) = (0.0);      (*value_1) = (0.1); */
    /* break; */
    /* case 23: /\* Sigmoidal_A0 *\//\* Sigmoidal Parameter in Sigma_0 and Sigma_1 (A_0) *\//\* 23 *\/ */
    /*   (*value_0) = (0.0);      (*value_1) = (0.1); */
    /* break; */
    /* case 24: /\* Sigmoidal_L0 *\//\* Sigmoidal Parameter in Sigma_0 and Sigma_1 (L_0) *\//\* 24 *\/ */
    /* (*value_0) = (1.0/100.0); (*value_1) = (10.0); */
    /* break; */
    /* case 25: /\* Sigmoidal_T0 *\//\* Sigmoidal Parameter in Sigma_0 and Sigma_1 (T_0) *\//\* 25 *\/ */
    /* (*value_0) = (2010.0); (*value_1) = (2015.0); */
    /* break; */
    default:
    printf("...(in Parameter_Space_Boundary_for_Trends.c) INVALID PARAMETER KEY (key = %d)\n", j);
    printf("...(in Parameter_Space_Boundary_for_Trends.c) VALID PARAMETER KEYS (from 0 to 20)\n");
    exit(0);
  }
}
