#include <MODEL.h>

#define EPSILON 1.0E-06

extern gsl_rng * r;

int M_O_D_E_L( Parameter_Table * Table )
{
  double x;

  int i,k, n;
  int I_Time, no_Patch;
  int Bad_Times;
  Parameter_Model * P;

  I_Time    = Table->T->I_Time;

  P = (Parameter_Model *)malloc( 1 * sizeof(Parameter_Model) );
  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (Table, P);
  Table->P  = P;
  printf(" Parameter_Model structure has been correctly allocated and initiated\n");

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  /* BEGIN : -------------------------------------------------------------------------
   * Definition Initial Condition (initializing 'Table->Vector_Model_Variables_Time_0' vector):
   */
  int MODEL_STATE_VARIABLES = K+1;
  Table->MODEL_STATE_VARIABLES = MODEL_STATE_VARIABLES;
  Table->Vector_Model_Variables_Time_0 = (double *)calloc( MODEL_STATE_VARIABLES, sizeof(double)); 

  if(Table->No_of_CELLS > 4) /* For instance, models DIFFUSION_AZTECA_4D and DIFFUSION_STOLLENBERG_4D */
    Initial_Condition_Centered_into_Parameter_Table (Table, Table->INITIAL_TOTAL_POPULATION);
    
  else if (Table->No_of_CELLS == 1)
    if(Table->TYPE_of_MODEL == 12 || Table->TYPE_of_MODEL == 13 || 
       Table->TYPE_of_MODEL == 14 || Table->TYPE_of_MODEL == 16   ) {
      Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
						   Table->TOTAL_No_of_FREE_CONSUMERS_TIME_0,
						   Table->TOTAL_No_of_HANDLING_CONSUMERS_TIME_0);
    }
    else {
      Initial_Condition_One_Single_Cell_into_Parameter_Table (Table,
							 Table->INITIAL_TOTAL_POPULATION,
							 Table->INITIAL_TOTAL_POPULATION);
    }
  else 
    Initial_Condition_All_Patches_the_Same_into_Parameter_Table (Table,
								 Table->INITIAL_TOTAL_POPULATION);
  /* END ----------------------------------------------------------------------------
   */

  /* BEGIN : -------------------------------------------------------------------------
   * Community Set Up
   */
  Community ** PATCH = (Community **)malloc( P->No_of_CELLS * sizeof(Community *) );
  Community_Allocation( PATCH, P ); 
  Community_Initialization (PATCH, P);
  /* The Parameter Model structure also keeps the three memmory addresses pointing to 
   * the Patch System, the Time Control structure, and the CPG structure to plot   
   */
  Table->Patch_System = PATCH;
  /* END ----------------------------------------------------------------------------
   */
			  							   
  Table->Vector_Model_Variables = (double *)calloc( MODEL_STATE_VARIABLES, sizeof(double) );

#if defined STATIONARY_POINT_REPRESENTATION  
  Table->Vector_Model_Variables_Stationarity = (double *)calloc( MODEL_STATE_VARIABLES,
								 sizeof(double) );
  Table->Vector_Model_Variables_MultiStability[0] = (double *)calloc( MODEL_STATE_VARIABLES,
								   sizeof(double) );
  Table->Vector_Model_Variables_MultiStability[1] = (double *)calloc( MODEL_STATE_VARIABLES,
								   sizeof(double) );
  Table->Vector_Model_Variables_MultiStability[2] = (double *)calloc( MODEL_STATE_VARIABLES,
								   sizeof(double) );

  #ifndef DIFFUSION_1R1C_2D
  /* B E G I N : Calculation of Stationary Points */  
  Fixed_Points_All( Table,                           
		    Table->Vector_Model_Variables_MultiStability[0],
		    Table->Vector_Model_Variables_MultiStability[1],
		    Table->Vector_Model_Variables_MultiStability[2],
		    EPSILON );

  /* The 'Function_to_Type_of_Stability()' function depends on the 
     JACOBIAN at the fixed point, which is only available for cetain models. 
     Please check the actual definition of the Jacobian matrix in the files 
     from the directory:  
     
        ./Definition_Numerical_Integration/ODE_Definitions/Include_ODE/

    of the type: 

        include.JAC_sys_[TYPE_of_MODEL].c
  */
  x = Function_to_Type_of_Stability_Double( Table );

  if (x == 0.0) printf("The Fixed Point is a Unstable Node\n");
  if (x == 1.0) printf("The Fixed Point is a Unstable Focus\n");
  if (x == 2.0) printf("The Fixed Point is a Stable Node\n");
  if (x == 3.0) printf("The Fixed Point is a Stable Focus\n");
  
  Print_Press_Key(1, 1, "Success: Fixed point and type of stability calculated!!!");  
  /*    E N D : --------------------------------- */
  #endif
#endif

  printf("\n");
  printf(" Entering deterministic dynamics. Parameter time dependencies will be\n");
  printf(" de-activated if -t4 0 (TYPE_of_TIME_DEPENDENCE = 0).\n");
  getchar();

  D_E_T_E_R_M_I_N_I_S_T_I_C___T_I_M_E___D_Y_N_A_M_I_C_S( Table ) ;

#if defined DIFFUSION_ECOEVO_PLANTS
  assert(Table->No_of_CELLS == 1);
  i = 0;
  double y_S, z_0, f_S;  
  
  y_S = Local_Population_Resources(i, Table->Vector_Model_Variables, Table);
  z_0 = (Table->K_R-y_S)/Table->K_R;
  
  for (k=0; k < Table->No_of_RESOURCES; k++) {
    n   = i*Table->LOCAL_STATE_VARIABLES + 2*k+1;
    
    y_S = Table->Vector_Model_Variables[n]/Table->K_R;
    
    f_S = Table->Beta_AP[k]/Table->Delta_AP[k] * z_0*Table->Eta_RP[k]*(1.0-2.0*Table->p_1) / (z_0 * Table->Eta_RP[k] + Table->Delta_RP[k]);

    printf("S=%d\t y=%.3f \t f=%.3f\t Beta=%.3f \t Eta=%.3f\n", 
            k+1,   y_S,      f_S, Table->Beta_AP[k], Table->Eta_RP[k]);
  } 

  Press_Key();
#endif 

#if defined CPGPLOT_REPRESENTATION
  /* Notice that j = TIMES now, as expected, since the program is just out
     from the loop:
        for( j = 1; j < TIMES; j++ ) { ... }
  */
  if(Table->TYPE_of_MODEL != 20) {   /* 20: MODEL = DIFFUSION_ECOEVO_PLANTS */
  ///    /* Temporal Evolution of the different output variables in separate subplots. 
  ///       For the ECOEVO_PLANTS model, only one dynamic bar plot  
  ///    */
    
    //  Parameter Table dependent costumized plotting is defined in ~/CPGPLOT/CPGPLOT_Parameter_Table/ files
    int TIMES           = Table->T->I_Time;
    int Input_Parameter = 0; /* The value of this model parameter appears in the title */
    // C_P_G___S_U_B___P_L_O_T_T_I_N_G ( Table, TIMES, Table->CPG->x_Time, Table->CPG->y_Time );
    C_P_G___S_U_B___P_L_O_T_T_I_N_G___C_U_S_T_O_M_I_Z_E_D___T_I_T_L_E ( Table,
   								                                                      TIMES,
   								                                                      Table->CPG->x_Time,
   								                                                      Table->CPG->y_Time,
   								                                                      Input_Parameter );
                        
  }
#endif

  free( Table->Vector_Model_Variables_Time_0);
  free( Table->Vector_Model_Variables );


#if defined STATIONARY_POINT_REPRESENTATION 
  // Fixed Points Calculations   
  #ifndef STO_REALIZATIONS
  /* De-allocation happens in MODEL_STO.c !!!
     If MODEL_STO.c is not used, then the following de-allocation should 
     be done here: 
  */
     free( Table->Vector_Model_Variables_MultiStability[0] );
     free( Table->Vector_Model_Variables_MultiStability[1] );
     free( Table->Vector_Model_Variables_MultiStability[2] );
     free( Table->Vector_Model_Variables_Stationarity );     
  #endif
#endif


  Community_Free(PATCH, P);
  free ( P );

  printf(" Deterministic dynamics successfully completed!!!\n");
  return(0);
}
