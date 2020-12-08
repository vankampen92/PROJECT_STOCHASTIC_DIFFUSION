 fprintf(fp, "\n");
 fprintf(fp, " P A R A M E T E R   S P A C E   D E F A U L T   V A L U E S :\n");

 fprintf(fp, " -sP %d\t No of parameters defining the actual dimension of the subparameter space within which the search will be performed: %d\n", No_of_PARAMETERS, No_of_PARAMETERS);

 fprintf(fp, "-I0 0: Parameter Key (%s) 0: \n", Table->Symbol_Parameters[0]);
 fprintf(fp, " [Parameter: 0]: -m0 %g -M0 %g -A0 %g -d0 %d\t [ min_Value[0] = %g, MAX_Value[0] = %g ]\t Accuracy[0] = %g, No of Points (along this axis) = %d\n",
	min_P0, MAX_P0, Acc_P0, d0, min_P0, MAX_P0, Acc_P0, d0);
 fprintf(fp, "-I1 1: Parameter Key (%s) 1: \n", Table->Symbol_Parameters[1]);
 fprintf(fp, " [Parameter: 1]: -m1 %g -M1 %g -A1 %g -d1 %d\t [ min_Value[1] = %g, MAX_Value[1] = %g ]\t Accuracy[1] = %g, No of Points (along this axis) = %d\n",
	min_P1, MAX_P1, Acc_P1, d1, min_P1, MAX_P1, Acc_P1, d1);
 fprintf(fp, "-I2 2: Parameter Key (%s) 2: \n", Table->Symbol_Parameters[2]);
 fprintf(fp, " [Parameter: 2]: -m2 %g -M2 %g -A2 %g -d2 %d\t [ min_Value[2] = %g, MAX_Value[2] = %g ]\t Accuracy[2] = %g, No of Points (along this axis) = %d\n",
	min_P2, MAX_P2, Acc_P2, d2, min_P2, MAX_P2, Acc_P2, d2);
 fprintf(fp, "-I3 3: Parameter Key (%s) 3: \n", Table->Symbol_Parameters[3]);
 fprintf(fp, " [Parameter: 3]: -m3 %g -M3 %g -A3 %g -d3 %d\t [ min_Value[3] = %g, MAX_Value[3] = %g ]\t Accuracy[3] = %g, No of Points (along this axis) = %d\n",
	min_P3, MAX_P3, Acc_P3, d3, min_P3, MAX_P3, Acc_P3, d3);
 fprintf(fp, "-I4 4: Parameter Key (%s) 4: \n", Table->Symbol_Parameters[4]);
 fprintf(fp, " [Parameter: 4]: -m4 %g -M4 %g -A4 %g -d4 %d\t [ min_Value[4] = %g, MAX_Value[4] = %g ]\t Accuracy[4] = %g, No of Points (along this axis) = %d\n",
  min_P4, MAX_P4, Acc_P4, d4, min_P4, MAX_P4, Acc_P4, d4);
 fprintf(fp, "-I5 5: Parameter Key (%s) 5: \n", Table->Symbol_Parameters[5]);
 fprintf(fp, " [Parameter: 5]: -m5 %g -M5 %g -A5 %g -d5 %d\t [ min_Value[5] = %g, MAX_Value[5] = %g ]\t Accuracy[5] = %g, No of Points (along this axis) = %d\n",
  min_P5, MAX_P5, Acc_P5, d5, min_P5, MAX_P5, Acc_P5, d5);
 fprintf(fp, "-I6 6: Parameter Key (%s) 6: \n", Table->Symbol_Parameters[6]);
 fprintf(fp, " [Parameter: 6]: -m6 %g -M6 %g -A6 %g -d6 %d\t [ min_Value[6] = %g, MAX_Value[6] = %g ]\t Accuracy[6] = %g, No of Points (along this axis) = %d\n",
 	min_P6, MAX_P6, Acc_P6, d6, min_P6, MAX_P6, Acc_P6, d6);
 fprintf(fp, " ... Parameter Space admits potentially more parameters\n"); 

// Parameter Space admits potentially more parameters. However, since Table->Symbol only allocates the
// the names of actual model parameters (from 0 to 6), the rest of this file should be commented out!
 /* fprintf(fp, "-I7 7: Parameter Key (%s) 7: Exponential Decay Parameter (f)\n", Table->Symbol[7]); */
 /* fprintf(fp, " [Parameter: 7]: -m7 %g -M7 %g -A7 %g -d7 %d\t [ min_Value[7] = %g, MAX_Value[7] = %g ]\t Accuracy[7] = %g, No of Points (along this axis) = %d\n", */
 /* 	min_P7, MAX_P7, Acc_P7, d7, min_P7, MAX_P7, Acc_P7, d7); */
 /* fprintf(fp, "-I8 8: Parameter Key (%s) 8: Sigmoidal Point (X_0)\n", Table->Symbol[8]); */
 /* fprintf(fp, " [Parameter: 8]: -m8 %g -M8 %g -A8 %g -d8 %d\t [ min_Value[8] = %g, MAX_Value[8] = %g ]\t Accuracy[8] = %g, No of Points (along this axis) = %d\n", */
 /* 	min_P8, MAX_P8, Acc_P8, d8, min_P8, MAX_P8, Acc_P8, d8); */
 /* fprintf(fp, "-I9 9: Parameter Key (%s) 9: Respiration (R_b)\n", Table->Symbol[9]); */
 /* fprintf(fp, " [Parameter: 9]: -m9 %g -M9 %g -A9 %g -d9 %d\t [ min_Value[9] = %g, MAX_Value[9] = %g ]\t Accuracy[9] = %g, No of Points (along this axis) = %d\n", */
 /* 	min_P9, MAX_P9, Acc_P9, d9, min_P9, MAX_P9, Acc_P9, d9); */

 /* fprintf(fp, "-I10 10: Parameter Key (%s) 10: Sea grass local birth rateColonization Rate\n", Table->Symbol[10]); */
 /* fprintf(fp, " [Parameter: 10]: -m10 %g -M10 %g -A10 %g -d10 %d\t [ min_Value[10] = %g, MAX_Value[10] = %g ]\t Accuracy[10] = %g, No of Points (along this axis) = %d\n", */
 /* min_P10, MAX_P10, Acc_P10, d10, min_P10, MAX_P10, Acc_P10, d10); */
 /* fprintf(fp, "-I11 11: Parameter Key (%s) 11: Sea grass death rate\n", Table->Symbol[11]); */
 /* fprintf(fp, " [Parameter: 11]: -m11 %g -M11 %g -A11 %g -d11 %d\t [ min_Value[11] = %g, MAX_Value[11] = %g ]\t Accuracy[11] = %g, No of Points (along this axis) = %d\n", */
 /* min_P11, MAX_P11, Acc_P11, d11, min_P11, MAX_P11, Acc_P11, d11); */
 /* fprintf(fp, "-I12 12: Parameter Key (%s) 12: Max No. of sea grass stems per unit area\n", Table->Symbol[12]); */
 /* fprintf(fp, " [Parameter: 12]: -m12 %g -M12 %g -A12 %g -d12 %d\t [ min_Value[12] = %g, MAX_Value[12] = %g ]\t Accuracy[12] = %g, No of Points (along this axis) = %d\n", */
 /* min_P12, MAX_P12, Acc_P12, d12, min_P12, MAX_P12, Acc_P12, d12); */
 /* fprintf(fp, "-I13 13: Parameter Key (%s) 13: Irradiance Level \n", Table->Symbol[13]); */
 /* fprintf(fp, " [Parameter: 13]: -m13 %g -M13 %g -A13 %g -d13 %d\t [ min_Value[13] = %g, MAX_Value[13] = %g ]\t Accuracy[13] = %g, No of Points (along this axis) = %d\n", */
 /* min_P13, MAX_P13, Acc_P13, d13, min_P13, MAX_P13, Acc_P13, d13); */
 /* fprintf(fp, "-I14 14: Parameter Key (%s) 14: Irradiance Half Constant\n", Table->Symbol[14]); */
 /* fprintf(fp, " [Parameter: 14]: -m14 %g -M14 %g -A14 %g -d14 %d\t [ min_Value[14] = %g, MAX_Value[14] = %g ]\t Accuracy[14] = %g, No of Points (along this axis) = %d\n", */
 /* min_P14, MAX_P14, Acc_P14, d14, min_P14, MAX_P14, Acc_P14, d14); */
 /* fprintf(fp, "-I15 15: Parameter Key (%s) 15: Exponential Decay Parameter (f)\n", Table->Symbol[15]); */
 /* fprintf(fp, " [Parameter: 15]: -m15 %g -M15 %g -A15 %g -d15 %d\t [ min_Value[15] = %g, MAX_Value[15] = %g ]\t Accuracy[15] = %g, No of Points (along this axis) = %d\n", */
 /* min_P15, MAX_P15, Acc_P15, d15, min_P15, MAX_P15, Acc_P15, d15); */
 /* fprintf(fp, "-I16 16: Parameter Key (%s) 16: Exponential Decay Parameter Beer-Lambert (gap)\n", Table->Symbol[16]); */
 /* fprintf(fp, " [Parameter: 16]: -m16 %g -M16 %g -A16 %g -d16 %d\t [ min_Value[16] = %g, MAX_Value[16] = %g ]\t Accuracy[16] = %g, No of Points (along this axis) = %d\n", */
 /* min_P16, MAX_P16, Acc_P16, d16, min_P16, MAX_P16, Acc_P16, d16); */
 /* fprintf(fp, "-I17 17: Parameter Key (%s) 17: Exponential Decay Parameter (f)\n", Table->Symbol[17]); */
 /* fprintf(fp, " [Parameter: 17]: -m17 %g -M17 %g -A17 %g -d17 %d\t [ min_Value[17] = %g, MAX_Value[17] = %g ]\t Accuracy[17] = %g, No of Points (along this axis) = %d\n", */
 /* min_P17, MAX_P17, Acc_P17, d17, min_P17, MAX_P17, Acc_P17, d17); */
 /* fprintf(fp, "-I18 18: Parameter Key (%s) 18: Sigmoidal Point (X_0)\n", Table->Symbol[18]); */
 /* fprintf(fp, " [Parameter: 18]: -m18 %g -M18 %g -A18 %g -d18 %d\t [ min_Value[18] = %g, MAX_Value[18] = %g ]\t Accuracy[18] = %g, No of Points (along this axis) = %d\n", */
 /* min_P18, MAX_P18, Acc_P18, d18, min_P18, MAX_P18, Acc_P18, d18); */
 /* fprintf(fp, "-I19 19: Parameter Key (%s) 19: Respiration (R_b)\n", Table->Symbol[19]); */
 /* fprintf(fp, " [Parameter: 19]: -m19 %g -M19 %g -A19 %g -d19 %d\t [ min_Value[19] = %g, MAX_Value[19] = %g ]\t Accuracy[19] = %g, No of Points (along this axis) = %d\n", */
 /* min_P19, MAX_P19, Acc_P19, d19, min_P19, MAX_P19, Acc_P19, d19); */

 fprintf(fp, " -sT %g\t TOLERANCE (minimization algorithm): %g\n", TOLERANCE, TOLERANCE);
 fprintf(fp, " -sN %d\t Max No of Iterations (minimization algorithm): %d\n", MAX_No_of_ITERATIONS,
       MAX_No_of_ITERATIONS );
