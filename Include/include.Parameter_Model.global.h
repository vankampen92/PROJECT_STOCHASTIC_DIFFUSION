int No_of_NEIGHBORS; 
int No_of_INDIVIDUALS;     
double Mu;                  /* Jumping Rate                    */
int No_of_CELLS;            /* No of CELLS                     */
int No_of_CELLS_X;          /* No of CELLS X Dimension         */
int No_of_CELLS_Y;          /* No of CELLS X Dimension         */
int No_of_RESOURCES;

double Lambda_R_0; 
double Delta_R_0; 
double Lambda_R_1; 
double Delta_R_1; 
int K_R;              /* Resource Carrying Capacity */ 

double Beta_R;        /* -H5 */ 
  
double Lambda_C_0;     /* -H5 */
double Delta_C_0;      /* -H6 */
double Lambda_C_1;     /* -H7 */
double Delta_C_1;      /* -H8 */ 
      
double Alpha_C_0;      /* -H9 */
double Nu_C_0;         /* -H10 */
  
double Chi_C_0;        /* -H11 */
double Eta_C_0;        /* -H12 */

double Mu_C;           /* -H13 */


int N_E;               /* -H14 */ /* Number of Energy Levels */
int f;                 /* -H15 */ /* Fecundity: Number of Offspring Individuals */
int i_0;               /* -H16 */ /* Energy Level at Maturity  */
double Beta_C;         /* -H17 */ /* Consummer Reproduction Rate */
int k_E;               /* -H18 */ /* 2* k_E is the resourse value in energy units */
double Theta_C;        /* -H19 */ /* Energy loss rate for maintenance */
double p_1;            /* -Hp1 */ /* Cooperation probability 1st position in the triplet */ 
double p_2;            /* -Hp2 */ /* Cooperation probability 2on position in the triplet */ 

/* Type of model */
#include <include.Type_of_Model.global.h>
