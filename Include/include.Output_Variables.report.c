fprintf(fp, "\n");
fprintf(fp, " Maximum number of output variables to plot: %d\n", OUTPUT_VARIABLES_MAXIMUM);
fprintf(fp, " In this run, the number of variable to plot is : ... ... > [PROGRAM] -n %d\n", SUB_OUTPUT_VARIABLES);
fprintf(fp, " To see the definition of the output variables, and so on,\n");
fprintf(fp, " see the C source file: ./definition_OutPut_Variables.c\n"); 
fprintf(fp, " The output variables that will be represented are: \n");
int i;

fprintf_Output_Variables(fp, P);

fprintf(fp, "\n ... ... ... ... ... ... > [PROGRAM] -n %d ", SUB_OUTPUT_VARIABLES); 
for (i = 0; i<SUB_OUTPUT_VARIABLES; i++){
  fprintf(fp, "-v%d %d ",  i, P->OUTPUT_VARIABLE_INDEX[i]); 
 }
fprintf(fp, "\n");





