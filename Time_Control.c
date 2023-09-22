#include <MODEL.h>

#include <include.Time_Control.extern.h>

void T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( Time_Control * Time, Parameter_Table * P, int I_Time)
{
  int i,j;

#if defined VERBOSE
  printf(" Time_Control is being allocated: \n");
  printf(" Several sets of %d output variables of length %d points will allocated\n",
	 P->SUB_OUTPUT_VARIABLES, I_Time);
#endif

  Time->AVE = (double **)calloc(P->SUB_OUTPUT_VARIABLES, sizeof(double *));
  for(i = 0; i<P->SUB_OUTPUT_VARIABLES; i++){
    Time->AVE[i]       = (double *)calloc( I_Time, sizeof(double) );
  }


  Time->VAR = (double **)calloc(P->SUB_OUTPUT_VARIABLES, sizeof(double *));
  for(i = 0; i<P->SUB_OUTPUT_VARIABLES; i++){
    Time->VAR[i]       = (double *)calloc( I_Time, sizeof(double) );
  }

  Time->time_DEF = (double *)calloc( I_Time, sizeof(double) );

  Time->count = (int *)calloc( I_Time, sizeof( int ) );

  Time->Time_Vector = (double *)calloc( I_Time, sizeof(double) );

  Time->summ = (double **)calloc(P->SUB_OUTPUT_VARIABLES, sizeof(double *));
  for(i = 0; i<P->SUB_OUTPUT_VARIABLES; i++){
    Time->summ[i]       = (double *)calloc( I_Time, sizeof(double) );
  }

#if defined VERBOSE
  printf(" Time_Control is being allocated: \n");
  printf(" %d output variables of length %d points will allocated\n",
	 P->SUB_OUTPUT_VARIABLES, I_Time);
#endif

  Time->summ_var   = (double **)calloc( P->SUB_OUTPUT_VARIABLES, sizeof(double *) );
  for (i=0; i<P->SUB_OUTPUT_VARIABLES; i++){
    Time->summ_var[i]   = (double *)calloc( I_Time, sizeof(double) );
  }

  Time->Rate = (Stochastic_Rate *)calloc( 1, sizeof(Stochastic_Rate) );

#if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
  Time->Vector_of_Rates = (double *)calloc(P->TOTAL_GRAND_No_of_EVENTS, sizeof(double)); 
#endif 

#if defined VERBOSE
  printf(" Time_Control is being allocated: \n");
  printf(" %d stochastic realizations of length %d points will allocated\n", P->Realizations, I_Time);
#endif

  Time->Accumulated_Variable = (double **)calloc( P->Realizations, sizeof(double *) );
  for (i=0; i<P->Realizations; i++){
    Time->Accumulated_Variable[i] = (double *)calloc( I_Time, sizeof(double) );
  }

  Time->Variable = (double ***)calloc(P->Realizations, sizeof(double **));
  for(i = 0; i<P->Realizations; i++){
    Time->Variable[i] = (double **)calloc( P->SUB_OUTPUT_VARIABLES, sizeof(double *) );
    for(j = 0; j<P->SUB_OUTPUT_VARIABLES; j++) {
      Time->Variable[i][j] = (double *)calloc( I_Time, sizeof(double ) );
    }
  }

  Time->Time_Vector_Real = (double **)calloc(P->Realizations, sizeof(double *));
  for(i = 0; i < P->Realizations; i++)
    Time->Time_Vector_Real[i] = (double *)calloc(I_Time, sizeof(double)); 
}

void  T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D ( Time_Control * Time, Parameter_Table * Table,
					                                    int I_Time )
{
  /* Setup for the vector of sampling times */
  int i;

  Time->EPSILON = EPSILON;
  Time->I_Time  = I_Time;

  Time->Time_0  = Time_0;
  Time->Time_1  = Time_1;
  
  Time->DISCARTING_EXTINCTIONS = 0; // DISCARTING_EXTINCTIONS to be included as input argument; 
  Time->Time_Scale_Unit = Time_Scale_Unit;
  Time->Delta_T = Delta_T; 
  Time->Rate->Stochastic_Time  = Time_0;
  Time->Realizations = Realizations;

  Time->TYPE_of_TIME_DEPENDENCE = TYPE_of_TIME_DEPENDENCE;
  
  for(i=0; i<I_Time; i++){
    Time->Time_Vector[i] = Time->Time_0 + (double)i * (Time->Time_1 - Time->Time_0)/(double)(I_Time-1);
    Time->count[i]       = 0;
  }
  /* so that Time_Vector[0] = Time_0 and  Time_Vector[I_Time-1] = Time_1. In this way, time series
     have I_Time points, where the first point always corresponds to Time_0 and the last to Time_1;
  */
  for(i = 0; i<Time->Realizations; i++) {
    Time->Time_Vector_Real[i][0] = Time->Time_0;
    Time->Time_Vector_Real[i][1] = Time->Time_0 + (double)(i+1) * (Time->Time_1 - Time->Time_0)/(double)(Time->Realizations);
  }

  Table->T = Time;
}

void T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( Time_Control * Time, Parameter_Table * P )
{
  int i,j;

  free (Time->time_DEF);

  for(i = 0; i<P->SUB_OUTPUT_VARIABLES; i++){
    free (Time->AVE[i]);
  }
  free (Time->AVE);
  for(i = 0; i<P->SUB_OUTPUT_VARIABLES; i++){
    free (Time->VAR[i]);
  }
  free (Time->VAR);

  free (Time->Time_Vector);

  free (Time->count );

  for(i = 0; i<P->SUB_OUTPUT_VARIABLES; i++){
    free (Time->summ[i]);
  }
  free (Time->summ);

  for (i=0; i<P->SUB_OUTPUT_VARIABLES; i++){
    free(Time->summ_var[i]);
  }
  free (Time->summ_var);

  for (i=0; i<P->Realizations; i++){
    free (Time->Accumulated_Variable[i]);
  }
  free (Time->Accumulated_Variable);

  free (Time->Rate);

#if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
  free(Time->Vector_of_Rates);  
#endif 

  for(i = 0; i<P->Realizations; i++){
    for(j = 0; j<P->SUB_OUTPUT_VARIABLES; j++) {
      free(Time->Variable[i][j]);
    }
    free(Time->Variable[i]);
  }
  free(Time->Variable);

  for(i = 0; i<P->Realizations; i++)
    free(Time->Time_Vector_Real[i]); 

  free(Time->Time_Vector_Real);
}

int Time_Control_AVE_VAR_SAVE_VARIABLES ( Parameter_Table * Table )
{
  /* This function is usually called after completing
     a set of stochastic realizations. It uses pervasively the Time_Control
     structure. As you see, it depends also on Parameter_Table structure.

     Notice that the corresponding vectors for each output variable summ[][]
     and summ_var[][] are used to calculate the averages and standard deviations
     respectively per time point.
  */
  int i,j,k, jj;
  int n        = Table->SUB_OUTPUT_VARIABLES;
  int I_TIMES  = Table->T->I_Time;
  double ave, var;
  Time_Control * Time = Table->T;

  FILE ** fp    = (FILE **)malloc( n * sizeof(FILE *) );
  char ** Files = (char **)calloc( n, sizeof(char *) );
  for(i=0; i < n; i++){
    jj = Table->OUTPUT_VARIABLE_INDEX[i];
    Files[i] = (char *)calloc( 50, sizeof(char) );
    Files[i][0]='\0';
    fitxer(Files[i], "output_VAR_", jj, ".dat");
  }

  /*Opening files */
  for(i=0; i < n; i++) fp[i]=fopen(Files[i],"w");

  /* Re-scaling times (daily), if necessary */
  /* ...   (using Time_Scale_Unit)          */
  /* Not done here                          */
  jj = 0;
  for (j=0; j<I_TIMES; j++){
    /* Saving the strain temporal evolution */
    if(Time->count[j] > 0){
      for(k=0; k < n; k++){
	ave = Time->summ[k][j]/(double)Time->count[j];
	var = Time->summ_var[k][j]/(double)Time->count[j] - ave*ave;
	if(var >= 0)
	  var = sqrt(var);
	else
	  var = 0.;
	fprintf(fp[k],"%g\t%g\t%g\n",Time->Time_Vector[j], ave, var);

	Time->AVE[k][jj] = ave;
	Time->VAR[k][jj] = var;
      }
      Time->time_DEF[jj] = Time->Time_Vector[j];
      jj++;
    }
  }
  /*Closing files and deallocating memmory */
  for(i=0; i < n; i++) { fclose(fp[i]); free (Files[i]); }
  free ( fp );
  free ( Files );

  return(jj);
}
