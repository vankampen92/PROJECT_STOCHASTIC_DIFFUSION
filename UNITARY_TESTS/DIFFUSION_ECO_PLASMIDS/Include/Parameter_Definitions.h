#define No_of_CELLS_SYSTEM          1        // Number of Local Communities in the Metacommunity 
#define INITIAL_CONDITION           1       // TYPE_of_INITIAL_CONDITION (1: Random; 0: Lower Bound) 

#define No_of_PLASMIDS_MAXIMUM      8     // 2     // 3      //  4     //   8     //   10 
#define No_of_PROFILES_MAXIMUM    256    //  4    //  8     //  16    //  256    //  1024   /* P = 2^(No_of_PLASMIDS_MAXIMUM)             */
#define No_of_STRAINS_MAXIMUM     100   //  10   // 100    //  100   //   100   //    100  /* S, number of Bacterial Types Maximum       */                                
#define No_of_RESOURCES_MAXIMUM 25600  //   40  //  800   //  1600  //  25600  //  102400 /* S * P */ /* P = 2^(No_of_PLASMIDS_MAXIMUM) */

/* No_of_RESOURCES_MAXIMUM is the Maximum No of Different strains IDs where S is the max number of 
   bacterial types or strain TYPEs */  

#define SYSTEM_SIZE                    1000    /* K_R        */ 
#define SPARSITY_PARAMETER              0.5    /* p_2        */          
#define CELL_DIVISION_RATE              5.0    /* Beta_R     */   
#define SEGREGATION_ERROR               0.1    /* p_1        */
#define BASAL_DEATH_RATE                1.0    /* Delta_R_0  */
#define STRESS_INDUCED_DEATH            0.5    /* Delta_R_1  */
#define COMPETITION_INDUCED_MORTALITY   0.5    /* Delta_C_0  */ 
#define COMMON_CONJUGATION_RATE         1.0    /* Lambda_R_1 */

#define DIFFUSION_RATE                  0.5    /* Mu         */

#define REPRODUCTION_COST               0.5    /* Alpha_C_0  */        
#define PLASMID_RESISTANCE              0.5    /* Nu_C_0     */
#define PLASMID_TRASMISSION_PROBABILITY 0.1    /* Chi_C_0    */

#define INITIAL_TIME                    0.0 
#define FINAL_TIME                    500.0
#define No_of_TIMES                   200
