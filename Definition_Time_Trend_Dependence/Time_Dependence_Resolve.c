#include <MODEL.h>

// #define SAVING_AN_INSTANCE_OF_THE_TREND

double Time_Dependence_Resolve(Parameter_Table * Table, int parameter, int pattern, double t)
{
  /* This sets the time-dependent parameter at its corresponding
  value at time t */
  int i,j,kk;

  double value;
  
  // printf(" Trend Control structure will be allocated: \n");
  Trend_Control * Tr = (Trend_Control *)malloc( 1 * sizeof(Trend_Control) );
  T_R_E_N_D___C_O_N_T_R_O_L___U_P_L_O_A_D( Tr, Table);
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

      default :
      printf(" Trend_Control.c: Type of Trend is not defined\n");
      printf(" Trend pattern can be only 0, 1, 2, 3 or 4.\n");
      printf(" but here Trend Type is %d\n", pattern); Press_Key();
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
  /* Time dependence is some parameters through "Sinusoidal Periodic Forcing":

        value = Mean_Value * ( 1.0 + Seasonal_Intensity * sin(2.0 * M_PI*t / Period) )

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
