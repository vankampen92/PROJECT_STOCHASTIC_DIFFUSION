typedef struct Parameter_Fittinginfo
{

  Observed_Data * Data; 
  Parameter_Table * Table; 
  Parameter_Space * Space;
  double (*Function)(const gsl_vector *, void *);
  int Verbose;
  int Minimization; 
  int Bounded_Parameter_Set;
  int TWO_PHASES;
  
}Parameter_Fitting;


