typedef struct Day_Month_RAIN_Info
{
  double day;
  double month;
  double rain;
}Day_Month_RAIN;

typedef struct Time_Dependence_Controlinfo
{
  /* This structure controls the way model parameters depend on time.
     Model parameters may depend on time through two different mechanisms:

     A. An "naturally driven" environmental covariate, which naturally fluctuates,
     incluencing model parameters. In this case, there is a time series that should be
     read from a file. There is also a 'response' function that translates environmental
     covariate fluctuation into model parameter fluctuations.

     B. Values from a file FILE_of_TIME_PARAMETERS

     C. An "artificially made" function of time, such as a sinusoidal, a linear trend,
     a triangular trend, a scholar-term seasonal forcing, or any other possible choice.
  */
  int No_of_TIMES;
  int No_of_EMPIRICAL_TIMES; 
  int TIME_DEPENDENT_PARAMETERS;
  int TYPE_of_TIME_DEPENDENCE;
  int TYPE_0_PARAMETERS;
  int TYPE_1_PARAMETERS;
  int TYPE_2_PARAMETERS;
  int No_of_COVARIATES;             /* Maximum number of co-variates */
  char ** Name_of_COVARIATES;

  double ** Dependent_Parameter;    /* Values for forced parameters at each time        */
  int * Index_Dependent_Parameters; /* Indeces of potentially forced model parameters */
  int * Forcing_Pattern_Parameters; /* Example: Index [] = {13, 3, 6, 15} */
                                    /* This means that only parameter 13, 3, 6, and 15 are forced */
                                    /* (by a forcing pattern of type 1).                 */
                                    /* Forcing_Patern [] = {1, 1, 1, 2}                  */
                                    /* This means parameter 13 is forced with pattern 1,  */
                                    /* while parameter 15 is forced with pattern 2.       */

  double ** COVARIATES;        /* Forcing parameters are, for instance, temperature or rainfall */
  char * FILE_of_COVARIATES;   /* This input file should have the following structure:          */
                               /* Time[0]                            Time[1]  ...               Time[No_of_TIMES-1] */
                               /* Value_CoVariate_1[0]  Value_CoVariate_1[1]  ...  Value_CoVariate_1[No_of_TIMES-1] */
                               /* Value_CoVariate_2[0]  Value_CoVariate_2[1]  ...  Value_CoVariate_2[No_of_TIMES-1] */

  char * FILE_of_TIME_PARAMETERS; /* This input file should have the following structure:          */
                               /* Time[0]                            Time[1]  ...               Time[No_of_TIMES-1] */
                               /* Value_Parameter_1[0]  Value_Parameter_1[1]  ...  Value_Parameter_1[No_of_TIMES-1] */
                               /* Value_Parameter_2[0]  Value_Parameter_2[1]  ...  Value_Parameter_2[No_of_TIMES-1] */
  double * Time_Vector;

  double ** TDC_Function_Params;  /* For each (potential) time dependent model parameter of type 2,
				                             this is an array with the parameters that control the 
				                             function that links time with true model parameters. 
				                          */
}Time_Dependence_Control;
