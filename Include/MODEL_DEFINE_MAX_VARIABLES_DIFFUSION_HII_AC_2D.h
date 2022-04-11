#define MODEL_PARAMETER_SPACE_MAXIMUM 5  /* Maximum Dimension for Parameter Space */

#define No_of_CELLS_MAXIMUM 40000         /* Let M be this maximum number of cells (i.e, 200 x 200) */

#define No_of_RESOURCES_MAXIMUM 2         /* S *//* Number of Local States */

#define MODEL_STATE_VARIABLES_MAXIMUM 80000  /* If M is the number of cells, then */
                                              /* Model dimension maximum is S * M */
#define OUTPUT_VARIABLES_TRUE_DERIVED 9       /* See definition_OutPut_Variables.c */

#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 11 /* Number Output Variables        */
                                         /* (other than MODEL_STATE_VARIABLES) */
                                         /* 11 = 2 (Local States) + 9 (True Derived) */
                                         /* See definition_OutPut_Variables.c  */

#define OUTPUT_VARIABLES_MAXIMUM 80011  /* MODEL_STATE_VARIABLES_MAXIMUM +
					   OUTPUT_VARIABLES_GENUINE_MAXIMUM */
