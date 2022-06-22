#include <MODEL.h>

extern gsl_rng * r;

extern int TYPE_of_TIME_DEPENDENCE;

int M_O_D_E_L___M_E( Parameter_Table * Table )
{
  int i,j,k, n;
  int I_Time, no_Patch;
  int Bad_Times;
  double t; 
  Time_Control * Time;
  Time_Dependence_Control * TDC; 

  Time = Table->T;
  TDC  = Table->TDC; 

  int No_of_CONFIGURATIONAL_STATES;
  int n_DIMENSION;
  int n_x, n_y, n_z;
  Model_Parameters_Master_Equation(Table,
				   &No_of_CONFIGURATIONAL_STATES,
				   &n_DIMENSION,
				   &n_x, &n_y, &n_z);
  
  Master_Equation * MEq = (Master_Equation *)calloc( 1, sizeof(Master_Equation) );
  /* MEq and Table will be two data structures that point to each other */
  Table->MEq = MEq;           
  MEq->Table = Table;
  
  Master_Equation_Allocation ( MEq,
			       No_of_CONFIGURATIONAL_STATES,
			       n_DIMENSION,
			       n_x, n_y, n_z );
 
  Master_Equation_Initialization ( MEq,
				   No_of_CONFIGURATIONAL_STATES,
				   n_DIMENSION,
				   n_x, n_y, n_z );
  
  Labels_for_Marginal_Probabilities( Table ); 
  
  /* Master Equation Numerical Integration                 */
  /* BEGIN: Core part (integration of the master equation) */
  printf("Entering Numerical Integration of the Master Equation...\n");   Press_Key();
  int ME_SYSTEM = master_equation_time_dynamics( Table );
  printf("Numerical Integration of the Master Equation succeeded!!!\n");  Press_Key();
  /*   END: ------------------------------------------ */

#if defined STOCHASTIC_REALIZATIONS  
 /* BEGIN : Representing distributions associated to output variables across 
            stochastic realization
 */  
#if defined CPGPLOT_REPRESENTATION
  int SAME_PLOT = 1;
  j = Table->T->I_Time - 1;
  assert(Table->MEq->n_DIMENSION <= 2);
  for ( i=0; i < Table->MEq->n_DIMENSION; i++ )
    C_P_G___E_M_P_I_R_I_C_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Table, j, i,
  							  Table->T->Time_Vector[j],
  							  SAME_PLOT );
  Press_Key();
#endif
 /*   END : Empirical Distribuiton Representation -----------------------*/
  
  /* Saving also the last time in a file */
  j = Table->T->I_Time - 1;
  assert(Table->MEq->n_DIMENSION <= 2);
  Saving_Empirical_Distribution_vs_ME_Numerical_Integration ( Table,
							      j, Table->T->Time_Vector[j] );
  
#endif
  
#if defined CPGPLOT_REPRESENTATION
  //  Parameter Table dependent costumized plotting is defined in
  //  ~/CPGPLOT/CPGPLOT_Parameter_Table/ files
  // Plotting the temporal evolution of the averages of the dynamical variables (which
  // should match the output variables for this work) calculated with the marginals from
  // the numerical integration of the master equation vs the densities calculated with
  // the ODE system.
  
  int TIMES           = Table->T->I_Time;
  //  Axes redefinition:
  Table->CPG->CPG_RANGE_X_0 = 0.0;  Table->CPG->CPG_RANGE_X_1 = Table->CPG->x_Time[TIMES-1];  
  Table->CPG->CPG_RANGE_Y_0 = 0.0;  Table->CPG->CPG_RANGE_Y_1 = (double)Table->No_of_INDIVIDUALS;
  Table->CPG->CPG_SCALE_Y   = 1;    Table->CPG->CPG_SCALE_X   = 1;
  SAME_PLOT = 0;
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T( SAME_PLOT,
						       Table, TIMES,
						       Table->CPG->x_Time,
						       Table->Matrix_Output_Variables);
  SAME_PLOT = 1;
  Press_Key();
  /* New colors and lines */
  Table->CPG->color_Index   = 3;
  Table->CPG->type_of_Width = 2;
  Table->CPG->type_of_Line  = 2;
  Table->CPG->type_of_Symbol = 5; 
  C_P_G___S_U_B___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T( SAME_PLOT,
						       Table, TIMES,
						       Table->CPG->x_Time,
						       Table->CPG->y_Time);
#endif

  FILE * fp = fopen("Data_ODE_vs_ME_Marginals_0.dat", "w");
  for(i=0; i<TIMES; i++)
    fprintf(fp, "%g\t%g\t%g\n", Table->CPG->x_Time[i],
	    Table->CPG->y_Time[0][i], Table->Matrix_Output_Variables[0][i]);  
  fclose(fp);
  
  if( Table->SUB_OUTPUT_VARIABLES > 1 ){
    fp = fopen("Data_ODE_vs_ME_Marginals_1.dat", "w");
    for(i=0; i<TIMES; i++)
      fprintf(fp, "%g\t%g\t%g\n", Table->CPG->x_Time[i],
	      Table->CPG->y_Time[1][i], Table->Matrix_Output_Variables[1][i]);  
    fclose(fp);
  }
  
  Master_Equation_Free ( MEq );

  return(0);
}

void C_P_G___E_M_P_I_R_I_C_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table,
							   int j,
							   int n, 
							   double Time_Current,
							   int SAME )
{
  /* This function plots marginal probability distributions at Time_Current, which 
     corresponds to the j-th time of the Time Vector stored in Table->T->Time_Vector. 

     As input, it takes Table. A member of Table is T, a pointer to a Time_Control structure 
     that stores all the necessary information to plot these distributions across stochastic 
     realizations. This distrubutions are empirical marginals of the output variables.  
     
     This function will plot these marginals, depending on an input parameter, an integer,
     n, defining the marginal to plot: 
     . n=0, first variable
     . n=1, second variable
     
     These empirical distributions can be compared with the marginal probability distributions 
     from the numerical integration of the master equation only if the first output variable
     is also the first dynamic variable and so on. If the output variables (in the same order)
     are NOT the same dynamic variables encoded in the master equation, this comparision will, 
     of course, fail. 
  */
  
  /* Notice that the probability distribution and marginals in the two dimensions are 
     associated to a Time_Current around the j-th time in Time_Vector[j] 
  */
  int Horizontal_Plot_Position, Vertical_Plot_Position;
  int i;
  int No_of_POINTS;
  static double Current_Time = 0.0;
  double Last_Time;

  Parameter_CPGPLOT * CPG = Table->CPG_STO; 
  Master_Equation * ME = Table->MEq;

  cpgslct(CPG->DEVICE_NUMBER);      /* Selecting Device */
  
  double * y;
  double * x;
  double * p; 
  
  if( ME->n_DIMENSION == 1 ) {
    
    y = (double *)calloc(Table->T->Realizations, sizeof(double));
    p = (double *)calloc(ME->n_x, sizeof(double));
    x = (double *)calloc(ME->n_x, sizeof(double));
    
    for(i=0; i<ME->n_x; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_x;
  }
  else if ( ME->n_DIMENSION == 2 ) {
    
    if ( n == 0 ) {
      y = (double *)calloc(Table->T->Realizations, sizeof(double));
      p = (double *)calloc(ME->n_x, sizeof(double));
      x = (double *)calloc(ME->n_x, sizeof(double));

      for(i=0; i<ME->n_x; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_x;
    }
    else if ( n == 1 ) { 
      y = (double *)calloc(Table->T->Realizations, sizeof(double));
      p = (double *)calloc(ME->n_y, sizeof(double));
      x = (double *)calloc(ME->n_y, sizeof(double));
      
      for(i=0; i<ME->n_y; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_y;
    }
    else {
      printf(" Error in master_equation_time_dynamics.c\n");
      printf(" The program will exit\n");
      exit(0);
    }
  }
  else {
    printf("Marginal probabilities for higher dimensions (n > %d) are not coded\n", n);
    printf("The program will exit\n");
    Press_Key();
    exit(0); 
  }

  int F; /* F: Number of errors at the j-th time across realizations; */ 
  int count, i_valid, valid_realizations; 
  count   = 0;
  i_valid = 0; 
  F       =  Table->T->Realizations - Table->T->count[j];
  for(i=0; i<Table->T->Realizations; i++) {
    
    if( Table->T->Variable[i][n][j] == 0.0 ) {
      if (count < F) count++;
      else {
	assert( Table->T->Variable[i][n][j] < No_of_POINTS );
	y[i_valid++] = Table->T->Variable[i][n][j];
      }
    }
    else {
      assert( Table->T->Variable[i][n][j] < No_of_POINTS );
      y[i_valid++] = Table->T->Variable[i][n][j];
    }
  }
  valid_realizations = i_valid;
  
  probability_distribution_from_stochastic_realizations( y, valid_realizations, 
							 p, No_of_POINTS );

  char * Plot_Title = (char *)calloc( 100, sizeof(char));
  char * Plot_Time  = (char *)calloc( 50, sizeof(char));
  char * Time_Eraser = (char *)calloc(50, sizeof(char));
  
  Last_Time     = Current_Time;
  Current_Time  = Time_Current; 
  sprintf(Plot_Time, "Time = %5.2f", Current_Time);
  sprintf(Time_Eraser, "Time = %5.2f", Last_Time);

  Plot_Title[0] ='\0';
  sprintf(Plot_Title, "Time = %5.2f", Current_Time);

  CPG->type_of_Line   = 3;
  CPG->type_of_Symbol = 5;
  CPG->type_of_Width  = 4;
  CPG->color_Index    = 3; 
  
  if (SAME > 0) {
    if (CPG->CPG__PANEL__X > 0 && CPG->CPG__PANEL__Y > 0 ){
      Horizontal_Plot_Position  = n%CPG->CPG__PANEL__X + 1;
      Vertical_Plot_Position    = n/CPG->CPG__PANEL__X + 1;
    }
    else {
      Vertical_Plot_Position    = n%abs(CPG->CPG__PANEL__Y) + 1;
      Horizontal_Plot_Position  = n/abs(CPG->CPG__PANEL__Y) + 1;;
    } 
    cpgpanl(Horizontal_Plot_Position, Vertical_Plot_Position);
    // printf("k = %d\t Horizontal Position = %d\t Vertical Position = %d\n",
    //	       k, Horizontal_Plot_Position, Vertical_Plot_Position);
    // Press_Key();
  }
  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME,
							No_of_POINTS, 
							x, p,
							ME->Marginal_Probability_Label[n], 
							"Probability", 
							Plot_Title,
							1, 1 );

  CPG->CPG_RANGE_X_0 = -0.5;  CPG->CPG_RANGE_X_1 = No_of_POINTS - 0.5;
  CPG->CPG_RANGE_Y_0 =  0.0;  CPG->CPG_RANGE_Y_1 = 1.0;  
  
  double x_Time_Position = 0.90 * (float)CPG->CPG_RANGE_X_1;
  double y_Time_Position = 0.90 * (float)CPG->CPG_RANGE_Y_1;

  cpgsch(2.0);
  cpgsci(0);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Time_Eraser);
  cpgsci(1);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Plot_Time);
  // cpgsch(char_Size);

  free(Plot_Title); free(Plot_Time); free(Time_Eraser);
  free(x); free(y); free(p); 
}

void Saving_Empirical_Distribution_vs_ME_Numerical_Integration ( Parameter_Table * Table,
								 int j,
								 double Time_Current)
{
  /* This function save the empirical probability distributions from stochastic 
     replicates vs the numerical integration of the Master Equaiton at Time_Current, which 
     corresponds to the j-th time of the Time Vector stored in Table->T->Time_Vector. 

     As input, it takes Table. A member of Table is T, a pointer to a Time_Control structure 
     that stores all the necessary information to plot these distributions across stochastic 
     realizations. 
     
     This empirical distributions can be compared with the probability distribution 
     from the numerical integration of the master equation only if the first output variable
     is also the first dynamic variable and so on. If the output variables are NOT the same 
     as dynamic variables encoded in the master equation (also in the same order), this 
     comparision will, of course, fail.

     For a 2D system, the output file will be a 4-columns file with this structure: 

     n  m  P_e  P_nm


     For a 1D system, the output file will be a 3-columns file with this structure: 

     n  P_e  P_n
     
     where P_e is the empirical probability and P_n is the result from performing 
     the numerical integration of the Master Equation until time Time_Current. 
  */
  
  /* Notice that the probability distribution is associated to a Time_Current 
     around the j-th time in Time_Vector[j] 
  */
  int i, n, m; 
  FILE * fp; 
  Master_Equation * ME = Table->MEq;

  double * y;
  double * x;
  double * P_n;
  double ** P_nm;
 
  int F; /* F: Number of errors at the j-th time across realizations; */ 
  int count, i_valid, valid_realizations; 

  fp = fopen("Probability_Distribution_Empirical_vs_Numerical.dat", "w");
  
  if( ME->n_DIMENSION == 1 ) {
   
    P_n = (double *)calloc(ME->n_x, sizeof(double));
    x   = (double *)calloc(Table->T->Realizations, sizeof(double));

    count   = 0;
    i_valid = 0; 
    F       =  Table->T->Realizations - Table->T->count[j];
    for(i=0; i<Table->T->Realizations; i++) {
      
      if( Table->T->Variable[i][0][j] == 0.0 ) {
	if (count < F) count++;
	else {
	  assert( Table->T->Variable[i][0][j] < ME->n_x );
	  x[i_valid++] = Table->T->Variable[i][0][j];
	}
      }
      else {
	assert( Table->T->Variable[i][0][j] < ME->n_x );
	x[i_valid++] = Table->T->Variable[i][0][j];
      }
    }
    valid_realizations = i_valid;
	    
    probability_distribution_from_stochastic_realizations( x, valid_realizations, 
							   P_n, ME->n_x );

    for(i = 0; i<ME->n_x; i++)
      fprintf(fp, "%d\t%g\t%g\n", i, P_n[i], ME->P_n[i]);

    free(x); free(P_n); 
  }
  
  else if ( ME->n_DIMENSION == 2 ) {

    P_nm = (double **)calloc(ME->n_x, sizeof(double *)); 
    for(i=0; i<ME->n_x; i++) 
      P_nm[i] = (double *)calloc(ME->n_y, sizeof(double));

    x   = (double *)calloc(Table->T->Realizations, sizeof(double));
    y   = (double *)calloc(Table->T->Realizations, sizeof(double));

    count   = 0;
    i_valid = 0; 
    F       =  Table->T->Realizations - Table->T->count[j];
    for(i=0; i<Table->T->Realizations; i++) {
      
      if( Table->T->Variable[i][0][j] == 0.0 && Table->T->Variable[i][1][j] == 0.0 )  {
	if (count < F)
	  count++;
	else {
	  x[i_valid] = Table->T->Variable[i][0][j];
	  y[i_valid] = Table->T->Variable[i][1][j];
	  i_valid++; 
	}
      }
      else {
	x[i_valid] = Table->T->Variable[i][0][j];
	y[i_valid] = Table->T->Variable[i][1][j];
	i_valid++;
      }
    }
    
    valid_realizations = i_valid;
    probability_distribution_from_stochastic_realizations_2D (x, y, valid_realizations,
							      P_nm, ME->n_x, ME->n_y);
    for(n = 0; n<ME->n_x; n++)
      for(m = 0; m<ME->n_y; m++) 
	fprintf(fp, "%d\t%d\t%g\t%g\n", n, m, P_nm[n][m], ME->P_nm[n][m]);

    free(x); free(y);
    
    for(i=0; i<ME->n_x; i++) 
      free(P_nm[i]);
    free(P_nm); 
  }
  
  else {
    printf("Marginal probabilities for higher dimensions (n > %d) are not coded\n", n);
    printf("The program will exit\n");
    Press_Key();
    exit(0); 
  }
  
  fclose(fp); 
}




