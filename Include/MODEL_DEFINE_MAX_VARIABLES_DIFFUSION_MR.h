#define MODEL_PARAMETER_SPACE_MAXIMUM 11  /* Maximum Dimension for Parameter Space */

#define No_of_CELLS_MAXIMUM 40000    /* Let M be this maximum number of cells (i.e, 200 x 200) */

#define No_of_RESOURCES_MAXIMUM 3    /* S *//* Number of local states */

#define MODEL_STATE_VARIABLES_MAXIMUM 120000  /* If M is the number of cells, then */
                                              /* Model dimension maximum is S * M */
#define OUTPUT_VARIABLES_TRUE_DERIVED 9       /* See definition_OutPut_Variables.c */

#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 12 /* Number Output Variables        */
                                         /* (other than MODEL_STATE_VARIABLES) */
                                         /* 12 = 3 (Local States Maximum) + 9 (True Derived) */
                                         /* See definition_OutPut_Variables.c  */

#define OUTPUT_VARIABLES_MAXIMUM 120012 /* MODEL_STATE_VARIABLES_MAXIMUM +
					   OUTPUT_VARIABLES_GENUINE_MAXIMUM */
