typedef struct Parameter_Fittinginfo
{

  Observed_Data * Data; 
  Parameter_Table * Table; 
  Parameter_Space * Space;
  /* It will be useful for this structure to point to 
     a function that depends on parameter values. This pointer will be 
     stored here 
  */ 
  double (*Function)(const gsl_vector *, void *);
  int Verbose;
  int Minimization; 
  int Bounded_Parameter_Set;
  int TWO_PHASES;
  int Iteration; 
  
}Parameter_Fitting;


