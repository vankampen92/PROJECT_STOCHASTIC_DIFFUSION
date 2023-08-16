#include <MODEL.h>

// #define SAVING_AN_INSTANCE_OF_THE_TREND

double Time_Dependence_Resolve(Parameter_Table * Table,
			       int parameter, int pattern, double t)
{
  /* Input Parameters: 
     . Table

     . parameter: index of the model parameter affected by time dependence 
     
     . n: order of the model parameter affected by time dependence in the vector 
     Table->TDC->Index_Dependent_Parameters (not included!!!)
     
     . pattern: type of trend that defines the fuction that links time and parameter
     values 

     . t: Current time
  */
  /* This sets the time-dependent parameter at its corresponding
  value at time t */
  int i,j,kk;
  
  double value;
  
  // printf(" Trend Control structure will be allocated: \n");
  Trend_Control * Tr = (Trend_Control *)malloc( 1 * sizeof(Trend_Control) );
  T_R_E_N_D___C_O_N_T_R_O_L___U_P_L_O_A_D( Tr, Table );
  Tr->TYPE_of_TREND  = pattern;
  // printf(" Trend_Control structure has been correctly allocated and initiated\n");

  switch (pattern)
  {
  case 0: //(for instance, sigmoidal)
    value = Time_Dependence_Sigmoidal_Function (Table, parameter, t);
    break;
    
  case 1: //(for instance, sinusoidal)
    value = Time_Dependence_Sinusoidal_Function (Table, parameter, t);
    break;
    
  case 2: //(first, linearly increasing and, then, decreasing)
    S_E_T_U_P___T_R_E_N_D___T_R_I_A_N_G_L_E ( Tr, Table->T, parameter );
    value = Triangular_Trend_Function (Table, t);
    break;
    
  case 3: //(linearly increasing)
    S_E_T_U_P___T_R_E_N_D___L_I_N_E_A_R ( Tr,  Table->T, parameter, +1);
    value = Linear_Trend_Function (Table, t);
    break;
    
  case 4: //(linearly decreasing)
    S_E_T_U_P___T_R_E_N_D___L_I_N_E_A_R ( Tr,  Table->T, parameter, -1);
    value = Linear_Trend_Function (Table, t);
    break;
    
  case 5: //(multi-step trend)
    S_E_T_U_P___T_R_E_N_D___M_U_L_T_I___S_T_E_P ( Tr, Table, parameter );
    value = Multi_Step_Function (Table, t);
    break;

  default :
    printf(" Trend_Control.c: Type of Trend is not defined\n");
    printf(" Trend pattern can be only 0, 1, 2, 3, 4 or 5\n");
    printf(" but here Trend Type is %d\n", pattern); Print_Press_Key(1,0,".");
    exit(0);
  }

  /* Saving 'value' at time t */
#if defined SAVING_AN_INSTANCE_OF_THE_TREND
  FILE *  fp;
  if ( (fp = fopen("Parameter_Trend.dat", "a")) == NULL ) {
    printf(" No file called Parameter_Trend.dat in current directorory\n");
    printf(" The program will exit\n");
    exit(0);
  }
  else {
    fprintf(fp, "%g\t%g\n", t, value);
  }
  fclose(fp);
#endif
  
  free(Tr);

  return(value);
}

double Time_Dependence_Sigmoidal_Function ( Parameter_Table * Table,
					    int parameter, double t)
{
  /* Time dependence is some parameters through "Sigmoidal Forcing":

       value = A0 /( 1.0 + exp(- Lambda * (t - T0)) );

     Obsolescent Function !!!
  */
  double value;
  double A0, T0, Lambda;

  value = 0.0;

  // Lambda = Table->Sigmoidal_L0;
  // A0     = Table->Sigmoidal_A0;
  // T0     = Table->Sigmoidal_T0;

  // value = A0 /( 1.0 + exp(- Lambda * (t - T0)) );

  return(value);
}

double Time_Dependence_Sinusoidal_Function ( Parameter_Table * Table,
					    int parameter, double t)
{
  /* Time dependence is some parameters through "Sinusoidal Periodic Forcing":

        value = Mean_Value * ( 1.0 + Seasonal_Intensity * sin(2.0 * M_PI*t / Period) )

	sinusoidal (of any possible kind of) forcing in some parameters

      Not implemented Function !!!
  */

  double value;
  double Seansonal_Intensity, Mean_Value, Variance_Value;

  value = 0.0;

  /*
     Seasonal_Intensity = Table->Sinusoidal_Seasonal_Intensity;
     Mean_Value         = Table->Sinusoidal_Mean_Value;
     Variance_Value     = Table->Sinusoidal_Variance_Value;
     Period             = Table->Sinusoidal_Period;

     value  = Sinusoidal_Oscillation(Time_Current,
                                     Seasonal_Intensity,
				     Mean_Value,
				     Variance_Value,
				     Period);
  */

  return(value);
}

double  Triangular_Trend_Function (Parameter_Table * P, double t)
{
  int Input_Parameter;
  double Slope;
  double value;
  double Time_0, Time_i, Time_1;

  Trend_Control * T = P->Tr;
  Time_0 = T->Tr_Time_0;
  Time_i = T->Tr_Time_i;
  Time_1 = T->Tr_Time_1;

  if ( t <= Time_0 ) value = T->Tr_value_0;
  else if ( t > Time_0 && t <= Time_i ) {
    Slope  = (T->Tr_value_i - T->Tr_value_0)/(Time_i - Time_0);
    value = T->Tr_value_0 + Slope * (t - Time_0);
  }
  else if ( t > Time_i && t < Time_1 ) {
    Slope  = (T->Tr_value_1 - T->Tr_value_i)/(Time_1 - Time_i);
    value = T->Tr_value_i + Slope * (t - Time_i);
  }
  else value = T->Tr_value_1;  /* if ( t >= Time_1 ) */

  return(value);
}

double  Linear_Trend_Function (Parameter_Table * P, double t)
{
  int Input_Parameter;
  double Slope;
  double value;
  double Time_0, Time_i, Time_1;

  Trend_Control * T = P->Tr;
  Time_0 = T->Tr_Time_0;
  Time_1 = T->Tr_Time_1;

  Slope  = (T->Tr_value_1 - T->Tr_value_0)/(Time_1 - Time_0);

  if ( t <= Time_0 )                    value = T->Tr_value_0;
  else if ( t > Time_0 && t <= Time_1 ) value = T->Tr_value_0 + Slope * (t - Time_0);
  else                                  value = T->Tr_value_1;

  return(value);
}

double Sinusoidal_Oscillation(double t, double E, double K_Mean, double K_Var, double T)
{
  double value;

  // E = gsl_ran_gaussian(r, K_var);
  /* This function returns a Gaussian random variate, with mean zero
     and standard deviation K_var */
  value = K_Mean *(1.0 + E * sin(2.*M_PI*t/T)); //gsl_ran_gaussian(r, K_Var);

  return (value);
}

double Multi_Step_Function (Parameter_Table * Table, double t)
{
  int i; 
  double value;
  double Central_Parameter_Value;
  
  int Input_Parameter = Table->Tr->Tr_Input_Parameter;
  int              N  = Table->Tr->Tr_No_of_Jumps;

  /* 
     This value is a reference value that has been introduced from command line, 
     which does not change during the whole execution of the code                
     This multistep function multiplies this reference value times a different 
     factor at different selected time intervals, which produces a discrete 
     step function 
  */

  /* Central parameter values does not change!!! They are initilized at the begining from 
     Table and copied into a Parameter_Model structure whose pointer is store in Table->P
  */
  Central_Parameter_Value = Parameter_Model_into_Vector_Entry ( Input_Parameter, Table->P );

  double * Values = (double *)calloc(2*N, sizeof(double)); 

  Values[0] = Table->Tr->Tr_value_0;
  Values[1] = Table->Tr->Tr_value_1;
  Values[2] = Table->Tr->Tr_value_2;
  Values[3] = Table->Tr->Tr_value_3;
  Values[4] = Table->Tr->Tr_value_4;
  Values[5] = Table->Tr->Tr_value_5;
  Values[6] = Table->Tr->Tr_value_6;
  Values[7] = Table->Tr->Tr_value_7;
  Values[8] = Table->Tr->Tr_value_8;
  Values[9] = Table->Tr->Tr_value_9;
  
  double * Jump_Times = (double *)calloc(2*N, sizeof(double)); 
  
  Jump_Times[0] = Table->Tr->Tr_Time_0;
  Jump_Times[1] = Table->Tr->Tr_Time_1;
  Jump_Times[2] = Table->Tr->Tr_Time_2;
  Jump_Times[3] = Table->Tr->Tr_Time_3;
  Jump_Times[4] = Table->Tr->Tr_Time_4;
  Jump_Times[5] = Table->Tr->Tr_Time_5;
  Jump_Times[6] = Table->Tr->Tr_Time_6;
  Jump_Times[7] = Table->Tr->Tr_Time_7;
  Jump_Times[8] = Table->Tr->Tr_Time_8;
  Jump_Times[9] = Table->Tr->Tr_Time_9;

  assert( Jump_Times[0] == Table->T->Time_0  && Jump_Times[N] == Table->T->Time_1 );
 
  for(i=0; i<N; i++) {
    if (t >= Jump_Times[i] && t < Jump_Times[i+1]) { 
      value = Values[i] * Central_Parameter_Value;
      break;
    }
  }

  free(Jump_Times); free(Values); 
  
  return(value);
}
