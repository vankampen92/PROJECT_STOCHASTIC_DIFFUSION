#define MODEL_PARAMETER_SPACE_MAXIMUM 21  /* Maximum Dimension for Parameter Space */

                                          /* DEPENDENT_PARAMETERS_MAXIMUM,
					     MODEL_PARAMETERS_MAXIMUM, and
					     MODEL_PARAMETER_SPACE_MAXIMUM
					     should always be defined to be
					     the same value
					  */

#define No_of_CELLS_MAXIMUM 2500          /* Let M be this maximum number of cells */

#define No_of_RESOURCES_MAXIMUM 121        /* S *//* Number of local states */
/* If the model includes energy levels, the S = (N_E + 1)*(N_E + 1);        */
/* where N_E is the number of energy levels                                 */

#define MODEL_STATE_VARIABLES_MAXIMUM 302500  /* If M is the number of cells, then */
                                              /* Model dimension maximum is S * M */
#define OUTPUT_VARIABLES_TRUE_DERIVED 9       /* See definition_OutPut_Variables.c */

#define OUTPUT_VARIABLES_GENUINE_MAXIMUM 130 /* Number Output Variables        */
                                         /* (other than MODEL_STATE_VARIABLES) */
                                         /* 130 = 121 (Local States Maximum) + 9 (True Derived) */
                                         /* See definition_OutPut_Variables.c  */

#define OUTPUT_VARIABLES_MAXIMUM 302630 /* MODEL_STATE_VARIABLES_MAXIMUM +
					   OUTPUT_VARIABLES_GENUINE_MAXIMUM */
