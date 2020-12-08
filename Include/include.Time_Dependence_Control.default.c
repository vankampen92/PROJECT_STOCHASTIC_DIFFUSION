TIME_DEPENDENT_PARAMETERS = DEPENDENT_PARAMETERS_MAXIMUM; 

No_of_COVARIATES   = 0; 
TYPE_0_PARAMETERS  = 0;
TYPE_1_PARAMETERS  = 4;
TYPE_2_PARAMETERS  = 0; 

dep_0 = 0;                
pat_0 = 0;

dep_1 = 1;                
pat_1 = 0;

dep_2 = 2;                
pat_2 = 0;

dep_3 = 3;                
pat_3 = 0;

dep_4 = 4;                
pat_4 = 0;

dep_5 = 5;                
pat_5 = 0;

dep_6 = 6;                
pat_6 = 0;

dep_7 = 7;
pat_7 = 0;

dep_9 = 9;
pat_9 = 0;

dep_10 = 10;
pat_10 = 0;          

dep_11 = 11;
pat_11 = 0;

dep_12 = 12;
pat_12 = 0;       

dep_13 = 13;
pat_13 = 0;

dep_14 = 14;
pat_14 = 0;

dep_15 = 15;
pat_15 = 0;

dep_16 = 16;
pat_16 = 0;

dep_17 = 17;
pat_17 = 0;

dep_18 = 18;
pat_18 = 0;

dep_19 = 19;
pat_19 = 0;

dep_20 = 20;
pat_20 = 0;          

dep_21 = 21;
pat_21 = 0;

dep_22 = 22;
pat_22 = 0;       

dep_23 = 23;
pat_23 = 0;

dep_24 = 24;
pat_24 = 0;

dep_25 = 25;
pat_25 = 0;

dep_26 = 26;
pat_26 = 0;

dep_27 = 27;
pat_27 = 0;

dep_28 = 28;
pat_28 = 0;

dep_29 = 29;
pat_29 = 0;

dep_30 = 30;
pat_30 = 0;          

dep_31 = 31;
pat_31 = 0;

dep_32 = 32;
pat_32 = 0;       

dep_33 = 33;
pat_33 = 0;

dep_34 = 34;
pat_34 = 0;

dep_35 = 35;
pat_35 = 0;

dep_36 = 36;
pat_36 = 0;

dep_37 = 17;
pat_37 = 0;

dep_38 = 18;
pat_38 = 0;

dep_39 = 19;
pat_39 = 0;


dependent_parameter = (int *)calloc(40, sizeof(int));
forcing_pattern     = (int *)calloc(40, sizeof(int));

// We have reserved 40 int space because this is what it is
// used in include.Time_Dependence_Control.default.aux.c !!!

#include <include.Time_Dependence_Control.default.aux.c>
