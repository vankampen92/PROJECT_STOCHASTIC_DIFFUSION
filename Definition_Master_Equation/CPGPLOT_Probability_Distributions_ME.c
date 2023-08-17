#include <MODEL.h>

/* CPGPLOTing Aux Functions: (Master_Equation_Functions.h)
    . void C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ()
    . void C_P_G___S_T_A_T_I_O_N_A_R_Y___D_I_S_T_R_I_B_U_T_I_O_N ()
    . void C_P_G___E_M_P_I_R_I_C_A_L___D_I_S_T_R_I_B_U_T_I_O_N ()
    . void C_P_G___T_H_E_O_R_E_T_I_C_A_L___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ()   
   SAVEing Aux Functions: 
    . void Saving_Empirical_Distribution_vs_ME_Numerical_Integration ()
    . void Saving_Marginal_Distribution()
*/

void C_P_G___T_H_E_O_R_E_T_I_C_A_L___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table,
							                                                                   int j,
							                                                                   int n, 
							                                                                   double Time_Current,
							                                                                   int SAME )
{
  /* This function plots time evolving theoretical marginal probability distributions. 
     As input, it takes Table. A member of Table is MEq, which stores all the necessary information 
     to plot the whole distribution and, of course, its theoretical marginals.

     Mathematical expressions for the theoretical marginals are not available for every model.
     So far, only available for DIFFUSION_HII_nD and only if all Nu's are the same. 

     Even if the dimension of the multivariate probibaility is low (n_DIMENSION is 1 or 2), 
     theoretical time evolving marginals are only correctly stored (if previously calculated) 
     in ME->MPD_T[k][n], where 'k' is the index of k-th time in Time_Vector[] and 'n' is 
     each of the dimension for every marginal. 

     Input Args: 
     . Table, 
     . j, index of the j-th time in Time_Vector[]
     . n, index of the marginal to plot
     . Time_Current
     . SAME, boolean argument telling whether reuse the same plot or not.   
     
     This function will plot the marginals, depending on an input parameter, an integer,
     n, defining the marginal to plot: 
     . n=0, first dimension
     . n=1, second dimension
     . n=2, third dimension, and so on.
  */
  
  /* Notice that the probability distribution and marginals in the two dimensions are 
     associated to a Time_Current around the j-th time in Time_Vector[j] 
  */
  int Horizontal_Plot_Position, Vertical_Plot_Position;
  int i;
  int No_of_POINTS;
  static double Current_Time = 0.0;
  double Last_Time;

  /* Checking if all Nu's are the same */
  assert_HLL_nD_Equal_Nus( Table );

  Parameter_CPGPLOT * CPG = Table->CPG_STO; 
  Master_Equation * ME = Table->MEq;

  cpgslct(CPG->DEVICE_NUMBER);      /* Selecting Device */
  
  double * y;
  double * x; 

  y = ME->MPD_T[j][n];
  x = (double *)calloc(ME->n_D[n], sizeof(double));   
  for(i=0; i<ME->n_D[n]; i++) 
    x[i] = (double)i;
    
  No_of_POINTS = ME->n_D[n]; 
  
  CPG->CPG_RANGE_X_0 = -0.5;  CPG->CPG_RANGE_X_1 = No_of_POINTS - 0.5;
  CPG->CPG_RANGE_Y_0 =  0.0;  CPG->CPG_RANGE_Y_1 = 1.0;  
  
  double x_Time_Position = 0.90 * (float)CPG->CPG_RANGE_X_1;
  double y_Time_Position = 0.90 * (float)CPG->CPG_RANGE_Y_1;

  char * Plot_Title = (char *)calloc( 100, sizeof(char));
  char * Plot_Time  = (char *)calloc( 50, sizeof(char));
  char * Time_Eraser = (char *)calloc(50, sizeof(char));
  
  Last_Time     = Current_Time;
  Current_Time  = Time_Current; 
  sprintf(Plot_Time, "Time = %5.2f", Current_Time);
  sprintf(Time_Eraser, "Time = %5.2f", Last_Time);

  Plot_Title[0] ='\0';
  sprintf(Plot_Title, "Time = %5.2f", Current_Time);

  CPG->type_of_Line   = 4;
  CPG->type_of_Symbol = 25;
  CPG->type_of_Width  = 8;
  CPG->color_Index    = 7; /* */

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
    // Print_Press_Key(1,0,".");
  }
  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME,
							                                          No_of_POINTS, 
							                                          x, y,
							                                          ME->Marginal_Probability_Label[n], 
							                                          "Probability", 
							                                          Plot_Title,
							                                          1, 1 );
  cpgsch(2.0);
  cpgsci(0);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Time_Eraser);
  cpgsci(1);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Plot_Time);
  // cpgsch(char_Size);

  free(Plot_Title); free(Plot_Time); free(Time_Eraser);
  free(x); 
}

void C_P_G___M_A_R_G_I_N_A_L___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table,
							                                           int j,
							                                           int n, 
							                                           double Time_Current,
							                                           int SAME )
{
  /* This function plots marginal probability distributions. As input, it takes Table. 
     A member of Table is MEq, which stores all the necessary information to plot the 
     whole distribution and, of course, its marginals.

     Input Args: 
     . Table, 
     . j, index of the j-th time in Time_Vector[]
     . n, index of the marginal to plot
     . Time_Current
     . SAME, boolean argument telling whether reuse the same plot or not.   
     
     This function will plot the marginals, depending on an input parameter, an integer,
     n, defining the marginal to plot: 
     . n=0, first dimension
     . n=1, second dimension
     . n=2, third dimension, and so on.
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
  
  if( ME->n_DIMENSION == 1 ) {
    assert(ME->n_DIMENSION == 1); 
    y = ME->P_n_Marginal;
    x = (double *)calloc(ME->n_x, sizeof(double));
    for(i=0; i<ME->n_x; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_x;
  }
  else if ( ME->n_DIMENSION == 2 ) {
    assert(ME->n_DIMENSION == 2);
    if ( n == 0 ) {
      y = ME->P_n_Marginal;
      x = (double *)calloc(ME->n_x, sizeof(double));
      for(i=0; i<ME->n_x; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_x;
    }
    else if ( n == 1 ) { 
      y = ME->P_m_Marginal;
      x = (double *)calloc(ME->n_y, sizeof(double));
      for(i=0; i<ME->n_y; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_y;
    }
    else {
      printf(" Error in file: CPGPLOT_Probability_Distributions_ME.c.\n");
      printf(" You are in two dimenison, but asking for the marginal probability\n");
      printf(" distribution for the higher dimension: %d\n", n);
      printf(" The program will exit\n");
      exit(0);
    }
  }
  else {
    y = ME->MPD[n];
    x = (double *)calloc(ME->n_D[n], sizeof(double));
    for(i=0; i<ME->n_D[n]; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_D[n]; 
  }
  
  CPG->CPG_RANGE_X_0 = -0.5;  CPG->CPG_RANGE_X_1 = No_of_POINTS - 0.5;
  CPG->CPG_RANGE_Y_0 =  0.0;  CPG->CPG_RANGE_Y_1 = 1.0;  
  
  double x_Time_Position = 0.90 * (float)CPG->CPG_RANGE_X_1;
  double y_Time_Position = 0.90 * (float)CPG->CPG_RANGE_Y_1;

  char * Plot_Title = (char *)calloc( 100, sizeof(char));
  char * Plot_Time  = (char *)calloc( 50, sizeof(char));
  char * Time_Eraser = (char *)calloc(50, sizeof(char));
  
  Last_Time     = Current_Time;
  Current_Time  = Time_Current; 
  sprintf(Plot_Time, "Time = %5.2f", Current_Time);
  sprintf(Time_Eraser, "Time = %5.2f", Last_Time);

  Plot_Title[0] ='\0';
  sprintf(Plot_Title, "Time = %5.2f", Current_Time);

  CPG->type_of_Line   = 4;
  CPG->type_of_Symbol = 25;
  CPG->type_of_Width  = 3;
  CPG->color_Index    = 6;

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
    // Print_Press_Key(1,0,".");
  }
  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME,
							                            No_of_POINTS, 
							                            x, y,
							                            ME->Marginal_Probability_Label[n], 
							                            "Probability", 
							                            Plot_Title,
							                            1, 1 );
  cpgsch(2.0);
  cpgsci(0);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Time_Eraser);
  cpgsci(1);
  cpgptxt(x_Time_Position, y_Time_Position, 0.0, 1.0, Plot_Time);
  // cpgsch(char_Size);

  free(Plot_Title); free(Plot_Time); free(Time_Eraser);
  free(x); 
}

void C_P_G___S_T_A_T_I_O_N_A_R_Y___D_I_S_T_R_I_B_U_T_I_O_N ( Parameter_Table * Table, 
							     int n, 
							     int SAME )
{
/* This function plots marginal theoretical probability distributions at stationarity.
     
  As input, it takes Table, the commanding structure of model parameters controling 
  model runs. SAME can be 0 or 1, depending on whether the plot have to be defined 
  from scratch (0) or not (1). 
     
  This function will plot these marginals, depending on an input parameter, an integer,
  n, defining the marginal to plot: 
     . n=0, first variable
     . n=1, second variable
     
  These theoretical distributions can be compared with the marginal probability distributions 
  from the numerical integration of the master equation when the final time is large enough, 
  and only if the first output variable is also the first dynamic variable and so on. If the 
  output variables (in the same order) are NOT the same dynamic variables encoded in the 
  master equation, this comparision will, of course, fail. 
*/
  
  int Horizontal_Plot_Position, Vertical_Plot_Position;
  int i;
  int No_of_POINTS;
  
  Parameter_CPGPLOT * CPG = Table->CPG_STO; 
  Master_Equation * ME = Table->MEq;

  cpgslct(CPG->DEVICE_NUMBER);      /* Selecting Device */
  
  double * x;
  double * p;                       /* Pointer to the marginal stationary probabilitiy 
				       distribution */
  
  if( ME->n_DIMENSION == 1 ) {
    
    p = ME->PS_n_Marginal; 
    x = (double *)calloc(ME->n_x, sizeof(double));
    
    for(i=0; i<ME->n_x; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_x;
  }
  else if ( ME->n_DIMENSION == 2 ) {
    
    if ( n == 0 ) {
      
      p = ME->PS_n_Marginal; 
      x = (double *)calloc(ME->n_x, sizeof(double));

      for(i=0; i<ME->n_x; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_x;
    }
    else if ( n == 1 ) { 
      
      p = ME->PS_m_Marginal; 
      x = (double *)calloc(ME->n_y, sizeof(double));
      
      for(i=0; i<ME->n_y; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_y;
    }
    else {
      printf(" Error in file: CPGPLOT_Probability_Distributions_ME.c.\n");
      printf(" You are in two dimensions, but asking for the marginal probability\n");
      printf(" distributjion the higher dimension: %d\n", n);
      printf(" The program will exit\n");
      exit(0);
    }
  }
  else {
    p = ME->MPD_S[n];
    x = (double *)calloc(ME->n_D[n], sizeof(double));
    for(i=0; i<ME->n_D[n]; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_D[n]; 
  }

  char * Plot_Title = (char *)calloc( 100, sizeof(char));
  Plot_Title[0] ='\0';
  // sprintf(Plot_Title, "Time = %5.2f", Current_Time);

  CPG->type_of_Line   = 3;
  CPG->type_of_Symbol = 5;
  CPG->type_of_Width  = 10;
  CPG->color_Index    = 5; 
  
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
    printf("n = %d\t Horizontal Position = %d\t Vertical Position = %d\n",
    	       n, Horizontal_Plot_Position, Vertical_Plot_Position);
    Print_Press_Key(1,0,".");
  }
  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME,
							                                          No_of_POINTS, 
							                                          x, p,
							                                          ME->Marginal_Probability_Label[n], 
							                                          "Probability", 
							                                          Plot_Title,
							                                          1, 1 );

  free(Plot_Title); 
  free(x); 
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
      printf(" Error in file: CPGPLOT_Probability_Distributions_ME.c.\n");
      printf(" You are in two dimensions, but asking for the marginal probability\n");
      printf(" distributjion the higher dimension: %d\n", n);
      printf(" The program will exit\n");
      exit(0);
    }
  }
  else {
      y = (double *)calloc(Table->T->Realizations, sizeof(double));
      p = (double *)calloc(ME->n_D[n], sizeof(double));
      x = (double *)calloc(ME->n_D[n], sizeof(double));
      
      for(i=0; i<ME->n_D[n]; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_D[n];
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
    // printf("n = %d\t Horizontal Position = %d\t Vertical Position = %d\n",
    //	       n, Horizontal_Plot_Position, Vertical_Plot_Position);
    // Print_Press_Key(1,0,".");
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

     The comparison between the full n-D distribution obtained through stochastic 
     replicates and the numerical integration of the n-D master equation is only 
     coded here for D < 3. For higher dimensions this function fails.
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
    printf(" Error in file: CPGPLOT_Probability_Distributions_ME.c.\n");
    printf(" You are asking for a comparison between the full n-D distribution obtained\n");
    printf(" through stochastic replicates and the numerical integration of the n-D master\n");
    printf(" equation for D = %d\n", ME->n_DIMENSION); 
    printf(" However, this function is only coded here for D < 3. For higher dimensions\n");
    printf(" this function fails.\n");
    printf(" The program will exit\n");
    exit(0); 
  }
  
  fclose(fp); 
}

void Saving_Marginal_Distribution(Parameter_Table * Table, int j, int n, double Time_Current)
{
  /* This function saves marginal probability distributions at a particular time, 
     Time_Current. As an input argument, it takes Table. A member of Table is MEq, 
     which stores all the necessary information to save the whole distribution and, 
     of course, its marginals.
     
     This function will plot the marginals, depending on an input parameter, an integer,
     n, defining the marginal to plot: 
     . n=0, first dimension
     . n=1, second dimension
     . n=2, third dimension (and so on, but ME->MPD[n] should have been calculated previously)
   */
  
  /* Notice that the probability distribution and marginals in the two dimensions are 
     associated to a Time_Current around the j-th time in Time_Vector[j] 
  */
  int i, k;
  int No_of_POINTS;
  
  Master_Equation * ME = Table->MEq;
  
  double * y;
  double * x; 
  
  if( ME->n_DIMENSION == 1 ) {
    assert(ME->n_DIMENSION == 1); 
    y = ME->P_n_Marginal;
    x = (double *)calloc(ME->n_x, sizeof(double));
    for(i=0; i<ME->n_x; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_x;
  }
  else if ( ME->n_DIMENSION == 2 ) {
    assert(ME->n_DIMENSION == 2);
    if ( n == 0 ) {
      y = ME->P_n_Marginal;
      x = (double *)calloc(ME->n_x, sizeof(double));
      for(i=0; i<ME->n_x; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_x;
    }
    else if ( n == 1 ) { 
      y = ME->P_m_Marginal;
      x = (double *)calloc(ME->n_y, sizeof(double));
      for(i=0; i<ME->n_y; i++) x[i] = (double)i;
      No_of_POINTS = ME->n_y;
    }
    else {
      printf(" Error in file: CPGPLOT_Probability_Distributions_ME.c.\n");
      printf(" You are in two dimensions, but asking for the marginal probability\n");
      printf(" distributjion the higher dimension: %d\n", n);
      printf(" The program will exit\n");
      exit(0);
    }
  }
  else {
    y = ME->MPD[n];
    x = (double *)calloc(ME->n_D[n], sizeof(double));
    for(i=0; i<ME->n_D[n]; i++) x[i] = (double)i;
    No_of_POINTS = ME->n_D[n]; 
  }

  char * NumDim   = (char *)calloc(10, sizeof(char));
  char * Marginal = (char *)calloc(100, sizeof(char));
  char * pFile;

  sprintf(NumDim, "%d", n);
  pFile = strcat(Marginal, "Marginal_Probability_");
  pFile = strcat(Marginal, NumDim);
  pFile = strcat(Marginal, "_Time_");

  printf("Saving Marginal Probability for Dimension %d at Time %g\n", n, Time_Current);  
  Saving_to_File_double(Marginal, x, y, No_of_POINTS, j);  
  
  free(Marginal); free(NumDim); 
  free(x); 
}
