typedef struct Parameter_Fittinginfo
{

  Observed_Data * Data; 
  Parameter_Table * Table; 
  Parameter_Space * Space;

  // Members to handle parametric configurations. 
  gsl_vector ** Solution;
  gsl_vector * Solution_Fitness;
  int * Solution_Order_Index;
  int No_of_SOLUTIONS; 
  int TOTAL_No_of_Fitting_Parameters;
  
  /* It will be useful for this structure to point to 
     a gsl function that depends on parameter values. This pointer will be 
     stored here 
  */ 
  double (*Function)(const gsl_vector *, void *);

  int Verbose;
  int Minimization; 
  int Bounded_Parameter_Set;
  int Ill_Defined_Function_Call; 
  int TWO_PHASES;
  int Iteration; 
  
}Parameter_Fitting;


