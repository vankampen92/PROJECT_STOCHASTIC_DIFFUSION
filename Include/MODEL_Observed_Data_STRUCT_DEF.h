typedef struct Observed_Datainfo
{
  /* A matrix of 'No_of_VARIABLES' observed variables, 'N', at a certain number of points, 
     'No_o_POINTS', for instance, time points, points along a spatial gradient, etc. 
  */
  double ** N;
  int No_of_VARIABLES;
  int No_of_POINTS; 
  char * Name;

  double * Time_Data_Vector;
  
}Observed_Data;
