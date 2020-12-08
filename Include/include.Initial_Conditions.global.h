/* Parameters Defining Initial Conditions for Human Population:
   see include.Initial_Conditions.default.c */
int i_POP, i_SP, n_MODEL_VARIABLE; 

int TYPE_of_INITIAL_CONDITION;

double INITIAL_TOTAL_POPULATION;

/* If this parameter is 1, then the initial total population is always rescaled to be:
 'INITIAT_TOTAL_POPULATION'. Otherwise, it is not!!!
*/
int RESCALING_INITIAL_TOTAL_POPULATION;

int No_of_IC;

int IC_0, IC_1, IC_2, IC_3, IC_4, IC_5, IC_6, IC_7, IC_8, IC_9;
int IC_d0, IC_d1, IC_d2, IC_d3, IC_d4, IC_d5, IC_d6, IC_d7, IC_d8, IC_d9;
double IC_MAX_0, IC_MAX_1, IC_MAX_2, IC_MAX_3, IC_MAX_4, IC_MAX_5, IC_MAX_6, IC_MAX_7, IC_MAX_8, IC_MAX_9;
double IC_min_0, IC_min_1, IC_min_2, IC_min_3, IC_min_4, IC_min_5, IC_min_6, IC_min_7, IC_min_8, IC_min_9;
double IC_Acc_0, IC_Acc_1, IC_Acc_2, IC_Acc_3, IC_Acc_4, IC_Acc_5, IC_Acc_6, IC_Acc_7, IC_Acc_8, IC_Acc_9;

int IC_10, IC_11, IC_12, IC_13, IC_14, IC_15, IC_16, IC_17, IC_18, IC_19;
int IC_d10, IC_d11, IC_d12, IC_d13, IC_d14, IC_d15, IC_d16, IC_d17, IC_d18, IC_d19;
double IC_MAX_10, IC_MAX_11, IC_MAX_12, IC_MAX_13, IC_MAX_14, IC_MAX_15, IC_MAX_16, IC_MAX_17, IC_MAX_18, IC_MAX_19;
double IC_min_10, IC_min_11, IC_min_12, IC_min_13, IC_min_14, IC_min_15, IC_min_16, IC_min_17, IC_min_18, IC_min_19;
double IC_Acc_10, IC_Acc_11, IC_Acc_12, IC_Acc_13, IC_Acc_14, IC_Acc_15, IC_Acc_16, IC_Acc_17, IC_Acc_18, IC_Acc_19;

int IC_20, IC_21, IC_22, IC_23, IC_24, IC_25, IC_26, IC_27, IC_28, IC_29;
int IC_d20, IC_d21, IC_d22, IC_d23, IC_d24, IC_d25, IC_d26, IC_d27, IC_d28, IC_d29;
double IC_MAX_20, IC_MAX_21, IC_MAX_22, IC_MAX_23, IC_MAX_24, IC_MAX_25, IC_MAX_26, IC_MAX_27, IC_MAX_28, IC_MAX_29;
double IC_min_20, IC_min_21, IC_min_22, IC_min_23, IC_min_24, IC_min_25, IC_min_26, IC_min_27, IC_min_28, IC_min_29;
double IC_Acc_20, IC_Acc_21, IC_Acc_22, IC_Acc_23, IC_Acc_24, IC_Acc_25, IC_Acc_26, IC_Acc_27, IC_Acc_28, IC_Acc_29;

int IC_30, IC_31, IC_32, IC_33, IC_34, IC_35, IC_36, IC_37, IC_38, IC_39;
int IC_d30, IC_d31, IC_d32, IC_d33, IC_d34, IC_d35, IC_d36, IC_d37, IC_d38, IC_d39;
double IC_MAX_30, IC_MAX_31, IC_MAX_32, IC_MAX_33, IC_MAX_34, IC_MAX_35, IC_MAX_36, IC_MAX_37, IC_MAX_38, IC_MAX_39;
double IC_min_30, IC_min_31, IC_min_32, IC_min_33, IC_min_34, IC_min_35, IC_min_36, IC_min_37, IC_min_38, IC_min_39;
double IC_Acc_30, IC_Acc_31, IC_Acc_32, IC_Acc_33, IC_Acc_34, IC_Acc_35, IC_Acc_36, IC_Acc_37, IC_Acc_38, IC_Acc_39;

int IC_40, IC_41, IC_42, IC_43;
int IC_d40, IC_d41, IC_d42, IC_d43;
double IC_MAX_40, IC_MAX_41, IC_MAX_42, IC_MAX_43;
double IC_min_40, IC_min_41, IC_min_42, IC_min_43;
double IC_Acc_40, IC_Acc_41, IC_Acc_42, IC_Acc_43;

double ** Ranges_IC;
double * Acc_IC;     //Accuracy for each parameter dimension
int    * d_IC;       //Discretization for each parameter dimension
int    * Index_IC;
