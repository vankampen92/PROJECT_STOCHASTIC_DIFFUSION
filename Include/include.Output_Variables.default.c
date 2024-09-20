SUB_OUTPUT_VARIABLES = OUTPUT_VARIABLES_MAXIMUM;

  variable_0 = 0;//2;
  variable_1 = 1; //16;
  variable_2 = 2; //4;
  variable_3 = 3; //5;
  variable_4 = 4; //6;
  variable_5 = 5; //8;

  variable_6  = 6;//2;
  variable_7  = 7; //16;
  variable_8  = 8; //4;
  variable_9  = 9; //5;
  variable_10 = 10; //6;
  variable_11 = 11; //8;

  variable_12 = 12;//2;

  variable_13 = 13; //16;
  variable_14 = 14; //4;
  variable_15 = 15; //5;
  variable_16 = 16; //6;
  variable_17 = 17; //8;

  variable_18 = 18;
  variable_19 = 19;
  variable_20 = 20;
  variable_21 = 21;
  variable_22 = 22;
  variable_23 = 23;
  variable_24 = 24;
  variable_25 = 25;
  variable_26 = 26;
  variable_27 = 27;
  variable_28 = 28;
  variable_29 = 29;
  variable_30 = 30;
  variable_31 = 31;
  variable_32 = 32;
  variable_33 = 33;
  variable_34 = 34;
  variable_35 = 35;
  variable_36 = 36;
  variable_37 = 37;
  variable_38 = 38;
  variable_39 = 39;

  variable_40 = 40;
  variable_41 = 41;
  variable_42 = 42;
  variable_43 = 43;
  variable_44 = 44;
  variable_45 = 45;
  variable_46 = 46;
  variable_47 = 47;
  variable_48 = 48;
  variable_49 = 49;

  variable_50 = 50;
  variable_51 = 51;
  variable_52 = 52;
  variable_53 = 53;
  variable_54 = 54;
  variable_55 = 55;
  variable_56 = 56;
  variable_57 = 57;
  variable_58 = 58;
  variable_59 = 59;

  Index_Output_Variables = (int *)calloc(OUTPUT_VARIABLES_MAXIMUM, sizeof(int));

// Alternatively, we may reserve less space, for instance, 
// 60 int space, and, then, make sure that this is what it is
// used in include.Output_Variables.default.aux.c !!!
// Index_Output_Variables = (int *)calloc(60, sizeof(int));

#include <include.Output_Variables.default.aux.c>
