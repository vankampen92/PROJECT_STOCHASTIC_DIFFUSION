 fprintf(fp, "\n");

 fprintf(fp, " \n Parameters Controling Trend Definition\n");

 fprintf(fp, " Trend Type ... ... ... ... ... ... ... ... ... ... ... . > PROGRAM -TY %d\n", TYPE_of_TREND);
 fprintf(fp, "   0. Trend Type = 0 TRIANGULAR TREND\n");
 fprintf(fp, "   1. Trend Type = 1 LINEAR TREND\n\n");
 
 fprintf(fp, " The initial starting time of the trend  ... ... . > PROGRAM -T0 %g\n", Tr_Time_0); 
 fprintf(fp, " The final time of the trend . ... ... .... .... . > PROGRAM -T1 %g\n", Tr_Time_1);
 fprintf(fp, " The intermediate value of time ... ... .. ... ... > PROGRAM -Ti %g\n\n", Tr_value_0);
 
 fprintf(fp, " The initial value of the parameter... ... ... . > PROGRAM -TI %g\n", Tr_value_0); 
 fprintf(fp, " The final value of the parameter ... ... .... . > PROGRAM -TF %g\n", Tr_value_1);
 fprintf(fp, " The intermediate value of parameter ... ... ... > PROGRAM -TM %g\n\n", Tr_value_i);
 
 
 fprintf(fp, " Trend Discretization (no of jumps) ... ... ... ... ... ..> PROGRAM -TJ %d\n", 
	Tr_No_of_Jumps);

 fprintf(fp, " Input model parameter key ... ... . ... ... ... ... ... ..> PROGRAM -TP %d\n", 
	Tr_Input_Parameter);
