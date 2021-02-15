No_of_INDIVIDUALS = 10000;     

No_of_CELLS = No_of_CELLS_MAXIMUM;

No_of_CELLS_X = 10; 

No_of_CELLS_Y = 10;

No_of_RESOURCES = No_of_RESOURCES_MAXIMUM; 

Mu      = 0.0;        /* -Hu  *//* Resources: Movement Per Capita Rate Between Patches */ 
Mu_C    = 0.0;        /* -H13 *//* Consumers: Movement Per Capita Rate Between Patches */ 

Lambda_R_0 = 0.0;      /* -H0 */
Delta_R_0  = 0.5;      /* -H1 */

K_R        = 5000.0;    /* -HK */ /* Resource Carrying Capacity */ 
Beta_R     = 1.5;      /* -H4 */ 
  
Lambda_C_0 = 10.0;     /* -H5 */
Delta_C_0  = 1.0;      /* -H6 */
      
Alpha_C_0 = 5.0;       /* -H9 */ /* Holling Type II Parameters */ 
Nu_C_0    = 1.0;      /* -H10 */
  
Chi_C_0 = 0.0;        /* -H11 */ /* Triplet formation  */
Eta_C_0 = 0.0;        /* -H12 */ /* Triplet degradation */

/* Additional 2nd Species */
Lambda_R_1 = 0.0;      /* -H2 */
Delta_R_1  = 0.0;      /* -H3 */

Lambda_C_1 = 0.0;      /* -H7 */
Delta_C_1 = 0.0;       /* -H8 */ 


#include <include.Type_of_Model.default.c>
