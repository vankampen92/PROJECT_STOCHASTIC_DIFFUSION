/* Initial Conditions parameters */

fprintf(fp, "\n");
fprintf(fp, " M O D E L   ( C O R E )   P A R A M E T E R S :\n");

#ifdef DIFFUSION_DRAG

fprintf(fp, " No of LOCAL POPULATIONS               ...> -HM   %d\n", No_of_CELLS);  
fprintf(fp, " No of LOCAL POPULATIONS (X dimension) ...> -HX   %d\n", No_of_CELLS_X);
fprintf(fp, " No of LOCAL POPULATIONS (Y dimension) ...> -HY   %d\n", No_of_CELLS_Y);
fprintf(fp, " Per Capita Movement Rate              ...> -Hu   %g\n", Mu);
fprintf(fp, " Total No of Indviduals                ...> -HN   %d\n", No_of_INDIVIDUALS);
fprintf(fp, " No of Species                         ...> -HS   %d\n", No_of_RESOURCES);

fprintf(fp, " External Immigration Rate (V)         ...> -H0   %g\n", Lambda_R_0);  
fprintf(fp, " Vegetation Decaying Rate  (V)         ...> -H1   %g\n", Delta_R_0);
fprintf(fp, " Vegetation Carrying Capicity (V)      ...> -HK   %d\n", K_R);
fprintf(fp, " Vegetation Per Capita Growth Rate (V) ...> -H4   %g\n", Beta_R);         /* -H4 */

fprintf(fp, " External Immigration Rate (R)         ---> -H5   %g\n", Lambda_C_0);     /* -H5 */
fprintf(fp, " Rat Per Capita Death Rate (R)         ---> -H6   %g\n", Delta_C_0);      /* -H6 */
fprintf(fp, " External Immigration Rate (G)         ---> -H7   %g\n", Lambda_C_1);     /* -H7 */
fprintf(fp, " Gull Per Capita Death Rate (G)        ---> -H8   %g\n", Delta_C_1);      /* -H8 */

fprintf(fp, "Rat Exploration/Attack Rate (R)       ---> -H9    %g\n", Alpha_C_0);      /* -H9 */
fprintf(fp, "Per Capita Movement Rate (R)          ---> -H13   %g\n", Mu_C);           /* -H13 */

fprintf(fp, "Gull Per Capita Reproduction Rate (G) ---> -H17   %g\n",  Beta_C);        /* -H17 */
fprintf(fp, "Energy Loss/Conversion Rate (R)       ---> -H19   %g\n",  Theta_C);       /* -H19 */
fprintf(fp, "Guano Conversion Factor (G)           ---> -Hp1   %g\n",  p_1);           /* -H20 */

#elif defined DIFFUSION_VRG

printf(fp, " No of LOCAL POPULATIONS               ...> -HM   %d\n", No_of_CELLS);
fprintf(fp, " No of LOCAL POPULATIONS (X dimension) ...> -HX   %d\n", No_of_CELLS_X);
fprintf(fp, " No of LOCAL POPULATIONS (Y dimension) ...> -HY   %d\n", No_of_CELLS_Y);
fprintf(fp, " Per Capita Movement Rate              ...> -Hu   %g\n", Mu);
fprintf(fp, " Total No of Indviduals                ...> -HN   %d\n", No_of_INDIVIDUALS);
fprintf(fp, " No of Species                         ...> -HS   %d\n", No_of_RESOURCES);

fprintf(fp, " External Immigration Rate (V)         ...> -H0   %g\n", Lambda_R_0);
fprintf(fp, " Vegetation Decaying Rate  (V)         ...> -H1   %g\n", Delta_R_0);
fprintf(fp, " Vegetation Carrying Capicity (V)      ...> -HK   %d\n", K_R);
fprintf(fp, " Vegetation Per Capita Growth Rate (V) ...> -H4   %g\n", Beta_R);         /* -H4 */

fprintf(fp, " External Immigration Rate (R)         ---> -H5   %g\n", Lambda_C_0);     /* -H5 */
fprintf(fp, " Rat Per Capita Death Rate (R)         ---> -H6   %g\n", Delta_C_0);      /* -H6 */
fprintf(fp, " External Immigration Rate (G)         ---> -H7   %g\n", Lambda_C_1);     /* -H7 */
fprintf(fp, " Gull Per Capita Death Rate (G)        ---> -H8   %g\n", Delta_C_1);      /* -H8 */

fprintf(fp, "Rat Exploration/Attack Rate (R)       ---> -H9    %g\n", Alpha_C_0);      /* -H9 */
fprintf(fp, "Per Capita Movement Rate (R)          ---> -H13   %g\n", Mu_C);           /* -H13 */

fprintf(fp, "Gull Per Capita Reproduction Rate (G) ---> -H17   %g\n",  Beta_C);        /* -H17 */
fprintf(fp, "Energy Loss/Conversion Rate (R)       ---> -H19   %g\n",  Theta_C);       /* -H19 */
fprintf(fp, "Guano Conversion Factor (G)           ---> -Hp1   %g\n",  p_1);           /* -H20 */

#else

fprintf(fp, " No of LOCAL POPULATIONS               ...> -HM   %d\n", No_of_CELLS);
fprintf(fp, " No of LOCAL POPULATIONS (X dimension) ...> -HX   %d\n", No_of_CELLS_X);
fprintf(fp, " No of LOCAL POPULATIONS (Y dimension) ...> -HY   %d\n", No_of_CELLS_Y);
fprintf(fp, " Per Capita Movement Rate (0)          ...> -Hu   %g\n", Mu);
fprintf(fp, " Total No of Indviduals                ...> -HN   %d\n", No_of_INDIVIDUALS);
fprintf(fp, " No of Species                         ...> -HS   %d\n", No_of_RESOURCES);

fprintf(fp, " External Immigration Rate (0)         ...> -H0   %g\n", Lambda_R_0);
fprintf(fp, " Resource Decaying Rate (0)            ...> -H1   %g\n", Delta_R_0);
fprintf(fp, " External Immigration Rate (1)         ...> -H2   %g\n", Lambda_R_1);
fprintf(fp, " Resource Decaying Rate (1)            ...> -H3   %g\n", Delta_R_1);
fprintf(fp, " Resource Patch Carrying Capicity      ...> -HK   %d\n", K_R);
fprintf(fp, " Resource Per Capita Growth Rate       ...> -H4   %g\n", Beta_R);         /* -H4 */

fprintf(fp, "Consumer External Immigration Rate (0) ---> -H5   %g\n", Lambda_C_0);     /* -H5 */
fprintf(fp, "Consumer Death Rate (0)                ---> -H6   %g\n", Delta_C_0);      /* -H6 */
fprintf(fp, "Consumer External Immigration Rate (0) ---> -H7   %g\n", Lambda_C_1);     /* -H7 */
fprintf(fp, "Consumer Death Rate (0)                ---> -H8   %g\n", Delta_C_1);      /* -H8 */

fprintf(fp, "Comsumer Attack Rate (0)              ---> -H9    %g\n", Alpha_C_0);      /* -H9 */
fprintf(fp, "Nu (One over the handling time)       ---> -H10   %g\n", Nu_C_0);         /* -H10 */
fprintf(fp, "Triplet formation rate (0)            ---> -H11   %g\n", Chi_C_0);        /* -H11 */
fprintf(fp, "Triplet desintegration rate (0)       ---> -H12   %g\n", Eta_C_0);        /* -H12 */

fprintf(fp, "Consumer movement rate (0)            ---> -H13   %g\n", Mu_C);           /* -H13 */

fprintf(fp, "Number of Energy Levels (0)           ---> -H14    %d\n", N_E);           /* -H14 */
fprintf(fp, "Fecundity: Per Capita Offspring No    ---> -H15   %d\n",  f);             /* -H15 */
fprintf(fp, "Energy Level at Maturity              ---> -H16   %d\n",  i_0);           /* -H16 */
fprintf(fp, "Consummer Reproduction Rate           ---> -H17   %g\n",   Beta_C);       /* -H17 */
fprintf(fp, "2* k_E: resourse value in energy units ---> -H18    %d\n", k_E);          /* -H18 */
fprintf(fp, "Maintenance Energy Loss Rate           ---> -H19   %g\n",  Theta_C);      /* -H19 */
fprintf(fp, "Cooperation probability 1st position   ---> -Hp1   %g\n", p_1);           /* -Hp1 */
fprintf(fp, "Cooperation probability 2on position   ---> -Hp2   %g\n", p_2);           /* -Hp2 */

fprintf(fp, "Propagule Establishment Rate           ---> -H20   %g\n", Eta_R);         /* -H20 */
#endif
