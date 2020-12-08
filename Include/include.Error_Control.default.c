MY_ERROR_HANDLER = 0; 

Err_0 = 0.1;
Err_1 = 0.1;
Err_2 = 0.1;
Err_3 = 0.1;
Err_4 = 0.1;
Err_5 = 0.1;
Err_6 = 0.1;
Err_7 = 0.1;
Err_8 = 0.1;
Err_9 = 0.1;
Err_10 = 0.1;
Err_11 = 0.1;
Err_12 = 0.1;
Err_13 = 0.1;
Err_14 = 0.1; 
Err_15 = 0.1; 
Err_16 = 0.1; 

TYPE_of_ERROR_FUNCTION = 0;

No_of_ERROR_PARAMETERS = 2;

Error_Parameter_1      = 0.1;

Error_Parameter_0      = 0.1;

E_Var_0 = 0;
E_d0 = 0;
E_MAX_0 = 1.0;
E_min_0 = 0.0;
E_Acc_0 = 0.01;

E_Var_1 = 1;
E_d1 = 0;
E_MAX_1 = 1.0;
E_min_1 = 0.0;
E_Acc_1 = 0.01;

E_Var_2 = 2;
E_d2 = 0;
E_MAX_2 = 1.0;
E_min_2 = 0.0;
E_Acc_2 = 0.01;

E_Var_3 = 3;
E_d3 = 0;
E_MAX_3 = 1.0;
E_min_3 = 0.0;
E_Acc_3 = 0.01;

E_Var_4 = 4;
E_d4 = 0;
E_MAX_4 = 1.0;
E_min_4 = 0.0;
E_Acc_4 = 0.01;

E_Var_5 = 5;
E_d5 = 0;
E_MAX_5 = 1.0;
E_min_5 = 0.0;
E_Acc_5 = 0.01;

E_Var_6 = 6;
E_d6 = 0;
E_MAX_6 = 1.0;
E_min_6 = 0.0;
E_Acc_6 = 0.01;

E_Var_7 = 7;
E_d7 = 0;
E_MAX_7 = 1.0;
E_min_7 = 0.0;
E_Acc_7 = 0.01;

E_Var_8 = 8;
E_d8 = 0;
E_MAX_8 = 1.0;
E_min_8 = 0.0;
E_Acc_8 = 0.01;

E_Var_9 = 9;
E_d9 = 0;
E_MAX_9 = 1.0;
E_min_9 = 0.0;
E_Acc_9 = 0.01;

E_Var_10 = 9;
E_d10 = 0;
E_MAX_10 = 1.0;
E_min_10 = 0.0;
E_Acc_10 = 0.01;

E_Var_11 = 9;
E_d11 = 0;
E_MAX_11 = 1.0;
E_min_11 = 0.0;
E_Acc_11 = 0.01;

E_Var_12 = 9;
E_d12 = 0;
E_MAX_12 = 1.0;
E_min_12 = 0.0;
E_Acc_12 = 0.01;

E_Var_13 = 9;
E_d13 = 0;
E_MAX_13 = 1.0;
E_min_13 = 0.0;
E_Acc_13 = 0.01;

E_Var_14 = 9;
E_d14 = 0;
E_MAX_14 = 1.0;
E_min_14 = 0.0;
E_Acc_14 = 0.01;

E_Var_15 = 9;
E_d15 = 0;
E_MAX_15 = 1.0;
E_min_15 = 0.0;
E_Acc_15 = 0.01;

E_Var_16 = 9;
E_d16 = 0;
E_MAX_16 = 1.0;
E_min_16 = 0.0;
E_Acc_16 = 0.01;

assert(OUTPUT_VARIABLES_MAXIMUM > 17); 

/* Ranges_Error = (double **)calloc(2, sizeof(double *)); */
/* Ranges_Error[0] = (double *)calloc(OUTPUT_VARIABLES_MAXIMUM, sizeof(double)); */
/* Ranges_Error[1] = (double *)calloc(OUTPUT_VARIABLES_MAXIMUM, sizeof(double)); */
/* Acc_Error       = (double *)calloc(OUTPUT_VARIABLES_MAXIMUM, sizeof(double)); */
/* d_Error         = (int *)calloc(OUTPUT_VARIABLES_MAXIMUM, sizeof(int)); */
/* Index_Error = (int *)calloc(OUTPUT_VARIABLES_MAXIMUM, sizeof(int)); */

Ranges_Error = (double **)calloc(2, sizeof(double *));
Ranges_Error[0] = (double *)calloc(17, sizeof(double));
Ranges_Error[1] = (double *)calloc(17, sizeof(double));
Acc_Error       = (double *)calloc(17, sizeof(double));
d_Error         = (int *)calloc(17, sizeof(int));
Index_Error = (int *)calloc(17, sizeof(int));

// We have reserved 17 int space because this is what it is
// used in include.Error_Control.default.aux.c !!!

#include <include.Error_Control.default.aux.c>
