#define MODEL_PARAMETER_SPACE_MAXIMUM 12 /* Maximum Dimension for Parameter Space */

#define No_of_CELLS_MAXIMUM 40000        /* Let M be the maximum number of cells (i.e, 200 x 200) */

#define No_of_RESOURCES_MAXIMUM 200      /* S * 2 */ /* Two states per species, RP and R */
                                         /* Maximum No of local states = S * 2, 
                                            where S is the max number of real species, resources or types 
                                            (S = 100), but each of them can be in two stages (RP or R) 
                                         */

#define MODEL_STATE_VARIABLES_MAXIMUM 8000000   /* If M is the number of cells, then */
                                                /* Model dimension maximum is 2 * S * M */
#define OUTPUT_VARIABLES_TRUE_DERIVED 9         /* See definition_OutPut_Variables.c */

#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 809 /* Number Output Variables        */
                                             /* (other than MODEL_STATE_VARIABLES) */
                                             /* 800 (Local States Maximum) + 9 (True Derived) */
                                             /* See definition_OutPut_Variables.c  */

#define OUTPUT_VARIABLES_MAXIMUM 8000809    /* MODEL_STATE_VARIABLES_MAXIMUM + OUTPUT_VARIABLES_GENUINE_MAXIMUM */
