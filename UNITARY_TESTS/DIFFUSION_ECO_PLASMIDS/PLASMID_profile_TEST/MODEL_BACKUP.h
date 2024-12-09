#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#ifndef GSL_HEADER_FILES
#define GSL_HEADER_FILES
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#endif

#define INITIAL_CONDITION         1       // TYPE_ofl_INITIAL_CONDITION: 0 (Lower bound); 1 (Random) 

#define No_of_PLASMIDS_MAXIMUM    2     //  2     // 3    //  8         
#define No_of_PROFILES_MAXIMUM    4    //   4    //  8   // 256     /* P = 2^(No_of_PLASMIDS_MAXIMUM)             */
#define No_of_STRAINS_MAXIMUM   100                                /* S, number of Bacterial Strains Maximum     */                                
#define No_of_RESOURCES_MAXIMUM 400  //   400  //  800 // 25600   /* S * P */ /* P = 2^(No_of_PLASMIDS_MAXIMUM) */
/* No_of_RESOURCES_MAXIMUM is the Maximum No of Different strains IDs 
   where S is the max number of bacterial types                   */  

#define No_of_CELLS_SYSTEM                1
#define SYSTEM_SIZE                    1000    /* K_R */ 
#define SPARSITY_PARAMETER              0.5    /* p_2 */          
#define CELL_DIVISION_RATE              1.0    /* Beta_R */   
#define SEGREGATION_ERROR               0.1    /* p_1 */
#define BASAL_DEATH_RATE                1.0    /* Delta_R_0 */
#define STRESS_INDUCED_DEATH            0.5    /* Delta_R_1 */
#define COMPETITION_INDUCED_MORTALITY   0.5    /* Delta_C_0 */ 
#define COMMON_CONJUGATION_RATE         1.0    /* Lambda_R_1 */
#define DIFFUSION_RATE                  0.5    /* Mu  */

#define REPRODUCTION_COST               0.5    /* Alpha_C_0 */        
#define PLASMID_RESISTANCE              0.5    /* Nu_C_0 */
#define PLASMID_TRASMISSION_PROBABILITY 0.1    /* Chi_C_0 */

#define INITIAL_TIME                    0.0 
#define FINAL_TIME                     10.0
#define No_of_TIMES                   100

#define VERBOSE
#define CURRENT_TIME_RDN_SEED

#if defined CPGPLOT_REPRESENTATION
/* This header file contains the function prototypes
 * of my hand-made cpgplot-based plotting procedures
 */
#include <CPGPLOT_HEADER.h>
#endif

#include <BP_StrainsProfiles.h>
