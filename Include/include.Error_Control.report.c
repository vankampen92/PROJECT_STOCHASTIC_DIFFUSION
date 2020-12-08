/* Initial Conditions parameters: Model Variables Initial Conditions */
fprintf(fp, "\n");
fprintf(fp, " E R R O R   C O N T R O L   A R G U M E N T S:\n");

fprintf(fp, " Label of 1st Output Variable affected by Error       ...> -eP0   %d\n", E_Var_0);  /*  8 */
fprintf(fp, " Label of 2on Output Variable affected by Error       ...> -eP1   %d\n", E_Var_1);  /*  9 */
fprintf(fp, " Label of 3rd Output Variable affected by Error       ...> -eP2   %d\n", E_Var_2);  /* 10 */
fprintf(fp, " Label of 4th Output Variable affected by Error       ...> -eP3   %d\n", E_Var_3);  /* 11 */

fprintf(fp, "  1st Output Variable Error Value      ...> -e0   %g\n", Err_0);  /*  8 */
fprintf(fp, "  2on Output Variable Error Value      ...> -e1   %g\n", Err_1);  /*  9 */
fprintf(fp, "  3rd Output Variable Error Value      ...> -e2   %g\n", Err_2);  /* 10 */
fprintf(fp, "  4th Output Variable Error Value      ...> -e3   %g\n", Err_3);  /* 11 */

fprintf(fp, " Type of Error Model: ...> -xn %d\n", TYPE_of_ERROR_FUNCTION);

fprintf(fp, " Number of error parametres affected by the minimization procedure: ...> -xN %d\n", No_of_ERROR_PARAMETERS);

fprintf(fp, "-eP0 0: Variable Index 0: \n");
fprintf(fp, "[Index of the 1st Variable (output variable): 0]: -em0 %g -eM0 %g\t [ min_Value[0] = %g, MAX_Value[0] = %g ]\t Accuracy[0] = %g, No of Points (along this axis) = %d\n",
	E_min_0, E_MAX_0, E_min_0, E_MAX_0, E_Acc_0, E_d0);
fprintf(fp, "-eP1 1: Variable Index 1: \n");
fprintf(fp, "[Index of the 2nd Variable (initial condition): 1]: -em0 %g -eM0 %g\t [ min_Value[0] = %g, MAX_Value[0] = %g ]\t Accuracy[0] = %g, No of Points (along this axis) = %d\n",
	E_min_1, E_MAX_1, E_min_1, E_MAX_1, E_Acc_1, E_d1);

fprintf(fp, "Note that accuracies and number of points along the axis can not be be changed from the command line, but modifying the corresponding default file!!!\n"); 

fprintf(fp, "Also, more error parameters can be included.. \n"); 

