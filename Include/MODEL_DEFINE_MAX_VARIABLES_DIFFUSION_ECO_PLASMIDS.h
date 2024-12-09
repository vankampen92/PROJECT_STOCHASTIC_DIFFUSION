#define MODEL_PARAMETER_SPACE_MAXIMUM 15 /* Maximum Dimension for Parameter Space */

// #define No_of_CELLS_MAXIMUM 2500         /* Let M be the maximum number of cells (i.e, 50 x 50) */
#define No_of_CELLS_MAXIMUM 100         /* Let M be the maximum number of cells (i.e, 10 x 10) */

// #define No_of_PLASMIDS_MAXIMUM 10        /* Number of Plasmids Maximum */ 
#define No_of_PLASMIDS_MAXIMUM 5        /* Number of Plasmids Maximum */

// #define No_of_PROFILES_MAXIMUM 1024      /* P = 2^(No_of_PLASMIDS_MAXIMUM)*/
#define No_of_PROFILES_MAXIMUM 32      /* P = 2^(No_of_PLASMIDS_MAXIMUM)*/

//#define No_of_STRAINS_MAXIMUM 100        /* S, number of Bacterial Strains Maximum           */
#define No_of_STRAINS_MAXIMUM 20        /* S, number of Bacterial Strains Maximum           */

// #define No_of_RESOURCES_MAXIMUM 102400   /* S * 2^10 */ /* 2^10 subpopulations per strain */
                                            /* Maximum No of local states = S * 2^No_of_PLASMIDS_MAXIMUM, 
                                               where S is the max number of strains (S = 100)  
                                            */
#define No_of_RESOURCES_MAXIMUM 640      /* S * 2^5 */ /* 2^10 subpopulations per strain */
                                         /* Maximum No of local states = S * 2^No_of_PLASMIDS_MAXIMUM, 
                                            where S is the max number of strains (S = 100)  
                                         */

// #define MODEL_STATE_VARIABLES_MAXIMUM 256000000  /* If M is the number of cells, then */
                                                    /* Model dimension maximum is S * P * M */
#define MODEL_STATE_VARIABLES_MAXIMUM 64000      /* If M is the number of cells, then */
                                                 /* Model dimension maximum is S * P * M */

#define OUTPUT_VARIABLES_TRUE_DERIVED 9          /* See definition_OutPut_Variables.c */

// #define OUTPUT_VARIABLES_GENUINE_MAXIMUM 102409 /* Number Output Variables        */
                                             /* (other than MODEL_STATE_VARIABLES) */
                                             /* 102400 (Local States Maximum) + 9 (True Derived) */
                                             /* See definition_OutPut_Variables.c  */
#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 649 /* Number Output Variables        */
                                               /* (other than MODEL_STATE_VARIABLES) */
                                               /* 102400 (Local States Maximum) + 9 (True Derived) */
                                               /* See definition_OutPut_Variables.c  */

// #define OUTPUT_VARIABLES_MAXIMUM 256102409    /* MODEL_STATE_VARIABLES_MAXIMUM + OUTPUT_VARIABLES_GENUINE_MAXIMUM */ 
#define OUTPUT_VARIABLES_MAXIMUM 64649        












/* MODEL_STATE_VARIABLES_MAXIMUM + OUTPUT_VARIABLES_GENUINE_MAXIMUM */ 




