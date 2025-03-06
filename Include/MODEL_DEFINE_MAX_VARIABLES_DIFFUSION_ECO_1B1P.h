#define MODEL_PARAMETER_SPACE_MAXIMUM 13/* Maximum Dimension for Parameter Space */

#define No_of_CELLS_MAXIMUM 10000         /* Let M be the maximum number of cells (i.e, 100 x 100) */

#define No_of_PLASMIDS_MAXIMUM 1        /* Number of Plasmids Maximum */

#define No_of_PROFILES_MAXIMUM 2      /* P = 2^(No_of_PLASMIDS_MAXIMUM)*/

#define No_of_STRAINS_MAXIMUM 1        /* S, number of Bacterial Strains Maximum           */

#define No_of_RESOURCES_MAXIMUM 2      /* Maximum No of local states = S * 2^No_of_PLASMIDS_MAXIMUM, 
                                            where S is the max number of strains (S = 100)  
                                       */

#define MODEL_STATE_VARIABLES_MAXIMUM 20000      /* If M is the number of cells, then */
                                                 /* Model dimension maximum is S * P * M */

#define OUTPUT_VARIABLES_TRUE_DERIVED 9          /* See definition_OutPut_Variables.c */

#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 11  /* Number Genuine Output Variables       
                                                (other than MODEL_STATE_VARIABLES) 
                                             */
                                             /* 11 = 2 (Local States Maximum) + 9 (True Derived) 
                                                See definition_OutPut_Variables.c  
                                              */
#define OUTPUT_VARIABLES_MAXIMUM 20011        




