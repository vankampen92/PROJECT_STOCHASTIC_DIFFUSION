#include <MODEL.h>

#if defined PARAMETRIC_CONFIGURATIONS
#if defined CPGPLOT_REPRESENTATION 

#include <include.CPG.extern.h>

void P_L_O_T_T_I_N_G___P_A_R_A_M_E_T_R_I_C___C_O_N_F_I_G_U_R_A_T_I_O_N_S (Parameter_Table * P, 
									  Parameter_Space * S, 
									  double ** X_0, 
									  int n_0, 
									  double ** X_1, 
									  int n_1)
{
  int i,j, k, n;
  int No_of_PANELS;
  
  double * x_Data = (double *)malloc( sizeof(double) * S->No_Parametric_Configurations );
  double * y_Data = (double *)malloc( sizeof(double) * S->No_Parametric_Configurations );

  assert( n_0 <= S->No_Parametric_Configurations );
  assert( n_1 <= S->No_Parametric_Configurations );

  No_of_PANELS = S->nPar * (S->nPar - 1) / 2; 
  //Parameter_CPGPLOT * CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( No_of_PANELS, 0, 0, "/CPS" );
  Parameter_CPGPLOT * CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( No_of_PANELS, 0, 0, "/XSERVE");
  printf(" Parameter_CPGPLOT plotting structure has been correctly allocated and initiated\n");

  FILE * f;
  char * p_Title;
  char * Y_label; 
  char * X_label; 
  char * Title     = (char *)malloc( sizeof(char) * 500 );
  char * File_Name = (char *)malloc( sizeof(char) * 500 );
  n = 0;
  for(i=0; i< S->nPar; i++){
    for(j=i+1; j< S->nPar; j++){
      
      k = S->I[i];  X_label = P->Name_Parameters[k];
      k = S->I[j];  Y_label = P->Name_Parameters[k];
      
      Title[0] = '\0';
      p_Title = strcat( Title, Y_label );    
      p_Title = strcat( Title, "_vs_" );
      p_Title = strcat( Title, X_label );    
      
      File_Name[0] = '\0';
      p_Title = strcat( File_Name, "Scattered_Plot_0_" );
      p_Title = strcat( File_Name, Title);
      p_Title = strcat( File_Name, ".dat" );
      f = fopen(File_Name, "w");
      printf("Name of file saved: %s\n", File_Name);

      for(k=0; k < n_0; k++) { 
	x_Data[k] = X_0[k][i];   
	y_Data[k] = X_0[k][j];
	fprintf(f, "%g\t%g\n", x_Data[k], y_Data[k]);
      }
      fclose(f);

      CPG->color_Index = 2; CPG->type_of_Symbol = 2;
      CPG->CPG_RANGE_X_0 = gsl_vector_get(S->P_min, i);
      CPG->CPG_RANGE_X_1 = gsl_vector_get(S->P_MAX, i);
      CPG->CPG_RANGE_Y_0 = gsl_vector_get(S->P_min, j);
      CPG->CPG_RANGE_Y_1 = gsl_vector_get(S->P_MAX, j);
      CPG->CPG_SCALE_X   = 1;
      CPG->CPG_SCALE_Y   = 1;
      CPGPLOT___X_Y___S_C_A_T_T_E_R_E_D___S_C_A_L_E ( CPG,
						      n_0, 
						      x_Data, y_Data,
						      X_label, Y_label, "",
						      CPG->CPG_SCALE_X, 
						      CPG->CPG_SCALE_Y );

      File_Name[0] = '\0';
      p_Title = strcat( File_Name, "Scattered_Plot_1_" );
      p_Title = strcat( File_Name, Title);
      p_Title = strcat( File_Name, ".dat" );
      f = fopen(File_Name, "w");
      for(k=0; k < n_1; k++) { 
	x_Data[k] = X_1[k][i];   
	y_Data[k] = X_1[k][j];
	fprintf(f, "%g\t%g\n", x_Data[k], y_Data[k]);
      }
      fclose(f);

      CPG->color_Index = 7; CPG->type_of_Symbol = 3; 
      CPGPLOT___X_Y___S_C_A_T_T_E_R_E_D___S_A_M_E___P_L_O_T ( CPG, 1,
							      n_1, 
							      x_Data, y_Data,
							      X_label, Y_label, "",
							      CPG->CPG_SCALE_X, 
							      CPG->CPG_SCALE_Y );

      if ( (++n)%(CPG->CPG__PANEL__X * CPG->CPG__PANEL__Y) == 0 ) Print_Press_Key(1,0,".");    
    }
  }
  
  free(Title); free(File_Name);
  free(x_Data); free(y_Data);
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( CPG, 0 );
  cpgclos();
 
  /* B E G I N : Saving a file with classified Parametric Configurations */
  f = fopen("Parametric_Configurations.dat", "w");
  for(j = 0; j < n_0; j++){
    fprintf(f, "%g\t", 0.0);
    for(i = 0; i < S->nPar; i++){
      fprintf(f, "%g\t", X_0[j][i]);
    }
    fprintf(f, "\n");
  }
  for(j = 0; j < n_1; j++){
    fprintf(f, "%g\t", 1.0);
    for(i = 0; i < S->nPar; i++){
      fprintf(f, "%g\t", X_1[j][i]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
  /*     E N D : ------------------------------------------------------ */
}
#endif
#endif
