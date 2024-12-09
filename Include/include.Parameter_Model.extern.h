extern int No_of_NEIGHBORS; 
extern int No_of_INDIVIDUALS;   /* Carrying Capacity for Juveniles */  
extern double Mu;               /* Jumping Rate                    */
extern int No_of_CELLS;         /* No of CELLS                     */
extern int No_of_CELLS_X;         /* No of CELLS X Dimension         */
extern int No_of_CELLS_Y;         /* No of CELLS Y Dimension         */
extern int No_of_RESOURCES;

extern int No_of_PLASMIDS; /* Used in DIFFUSION_ECO_PLASMIDS */
extern int No_of_STRAINS;  /* Used in DIFFUSION_ECO_PLASMIDS */

extern double Lambda_R_0;    /* -H0 */
extern double Delta_R_0;     /* -H1 */
extern double Lambda_R_1;    /* -H2 */
extern double Delta_R_1;     /* -H3 */
extern int K_R;              /* -HK: Resource Carrying Capacity */ 

extern double Beta_R;        /* -H5 */ 
  
extern double Lambda_C_0;     /* -H5 */
extern double Delta_C_0;      /* -H6 */
extern double Lambda_C_1;     /* -H7 */
extern double Delta_C_1;      /* -H8 */ 
      
extern double Alpha_C_0;      /* -H9 */
extern double Nu_C_0;         /* -H10 */
  
extern double Chi_C_0;        /* -H11 */
extern double Eta_C_0;        /* -H12 */

extern double Mu_C;           /* -H13 */

extern int N_E;               /* -H14 */ /* Number of Energy Levels */
extern int f;                 /* -H15 */ /* Fecundity: Number of Offspring Individuals */
extern int i_0;               /* -H16 */ /* Energy Level at Maturity  */
extern double Beta_C;         /* -H17 */ /* Consummer Reproduction Rate */
extern int k_E;               /* -H18 */ /* 2* k_E is the resourse value in energy units */
extern double Theta_C;        /* -H19 */ /* Energy loss rate for maintenance */
extern double p_1;     /* -Hp1 */ /* Cooperation probability 1st position in the triplet */ 
extern double p_2;     /* -Hp2 */ /* Cooperation probability 2on position in the triplet */ 

extern double Eta_R;           /* -H20 */ /* Propagule establishment rate */

#include <include.Type_of_Model.extern.h>
