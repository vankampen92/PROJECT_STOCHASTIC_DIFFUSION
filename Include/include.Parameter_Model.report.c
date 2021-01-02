/* Initial Conditions parameters */

fprintf(fp, "\n");
fprintf(fp, " M O D E L   ( C O R E )   P A R A M E T E R S :\n");
fprintf(fp, " No of LOCAL POPULATIONS               ...> -HM   %d\n", No_of_CELLS);  
fprintf(fp, " No of LOCAL POPULATIONS (X dimension) ...> -HX   %d\n", No_of_CELLS_X);
fprintf(fp, " No of LOCAL POPULATIONS (Y dimension) ...> -HY   %d\n", No_of_CELLS_Y);
fprintf(fp, " Per Capita Movement Rate              ...> -Hu   %g\n", Mu);
fprintf(fp, " Total No of Indviduals                ...> -HN   %d\n", No_of_INDIVIDUALS);
fprintf(fp, " No of Species                         ...> -HS   %d\n", No_of_RESOURCES);

fprintf(fp, " External Immigration Rate (0)         ...> -H0   %g\n", Lambda_R_0);  
fprintf(fp, " Resource Decaying Rate (0)            ...> -H1   %g\n", Delta_R_0);
fprintf(fp, " External Immigration Rate (1)         ...> -H2   %g\n", Lambda_R_1);
fprintf(fp, " Resource Decaying Rate (1)            ...> -H3   %g\n", Delta_R_1);
fprintf(fp, " Resource Patch Carrying Capicity      ...> -HK   %d\n", K_R);


