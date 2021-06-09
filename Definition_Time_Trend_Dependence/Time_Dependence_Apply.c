#include <MODEL.h>

void Time_Dependence_Apply(Parameter_Table * Table, double t)
{
  /* This sets the time-dependent parameter at its corresponding
  value at time t */
  int i, j, kk, N;
  double value;
  
  N  = Table->TDC->TIME_DEPENDENT_PARAMETERS;
  
  j  = Integer_Position_of_a_Time(Table, t);
  assert( j >= 0 && j < Table->TDC->No_of_TIMES);
  
  for(i=0; i<N; i++) {
      kk = Table->TDC->Index_Dependent_Parameters[i];
    
      assert( kk >= 0 && kk < MODEL_PARAMETERS_MAXIMUM);
      
      value = Table->TDC->Dependent_Parameter[i][j];

      AssignVectorEntry_to_Structure(Table, kk, value);
  }
}

void Time_Dependence_Apply_Optimized(Parameter_Table * Table, int j)
{
  /* This sets the time-dependent parameter at its corresponding
  value at time t */
  int i, kk, N;
  double value;
  
  N  = Table->TDC->TIME_DEPENDENT_PARAMETERS;
  
  assert( j >= 0 && j < Table->TDC->No_of_TIMES);

  assert(Table->T->I_Time == Table->TDC->No_of_TIMES); 
  
  for(i=0; i<N; i++) {
      kk = Table->TDC->Index_Dependent_Parameters[i];
    
      assert( kk >= 0 && kk < MODEL_PARAMETERS_MAXIMUM);
      
      value = Table->TDC->Dependent_Parameter[i][j];

      AssignVectorEntry_to_Structure(Table, kk, value);
  }
}

int Integer_Position_of_a_Time(Parameter_Table * Table, double t)
{
  int n, position;
  int position_found;
  double t_0, t_1;

  t_0 = Table->TDC->Time_Vector[0];
  t_1 = Table->TDC->Time_Vector[1];
  position_found = 0; n = 0;
  while (position_found == 0){
    if( t >= t_0 && t < (t_0 + (t_1-t_0)/2.0) ) {
      position_found = 1;
      position = n;
    }
    else if ( t >= (t_0 + (t_1-t_0)/2.0) && t <= t_1 ) {
      position_found = 1;
      position = n+1;
    }
    else {
      n++;

      assert( n < (Table->TDC->No_of_TIMES-1) );

      t_0 = Table->TDC->Time_Vector[n];
      t_1 = Table->TDC->Time_Vector[n+1];
    }
  }
  return(position);
}
