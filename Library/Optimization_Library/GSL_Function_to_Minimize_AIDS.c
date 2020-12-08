#include <MODEL.h>

int M_O_D_E_L( Parameter_Table * );

double GSL_Function_to_Minimize_AIDS( const gsl_vector * x, void * Par )
{
  /* This GSL function allows an evaluation of function to minimize

     Output value:
     . Value
  */
  static double Value; 
  double Value_0, Value_1, Value_2;
  int j;

  Parameter_Fitting * F = (Parameter_Fitting *)Par;

  int No_of_POINTS         = F->Data->No_of_POINTS;
  int No_of_VARIABLES      = F->Data->No_of_VARIABLES;
  double ** Data           = F->Data->N;
  double ** Theory         = F->Table->Matrix_Output_Variables;
  Parameter_Table * Table  = F->Table; 
  Parameter_Space * Space  = F->Space;
  int No_of_PARAMETERS     = F->Space->No_of_PARAMETERS;
  int x_is_BOUNDED; 
  
  x_is_BOUNDED = 1; // Checking_for_Parameter_Boundaries( F, x );

  if( x_is_BOUNDED == 1 ) { 
						    
    Vector_Entries_into_Parameter_Table ( x, Table,
					  Space->Parameter_Index, No_of_PARAMETERS );
    
    if(Table->No_of_IC > 0) 
      Vector_Entries_into_Parameter_Table_Initial_Condition ( x, Table, 
							      Table->IC_Space->Parameter_Index,
							      No_of_PARAMETERS,
							      Table->IC_Space->No_of_PARAMETERS);
    int State = M_O_D_E_L(Table);
    
    assert(State == GSL_SUCCESS);   
    
    assert(No_of_VARIABLES == 3);  
    
    Value_0 = 0.0; 
    for( j=0; j<No_of_POINTS; j++ ) 
      if (Data[0][j] > 0.0) 
	Value_0 += (Data[0][j] - Theory[0][j]) * (Data[0][j] - Theory[0][j]);
    
    Value_1 = 0.0; 
    for( j=0; j<No_of_POINTS; j++ ) 
      if (Data[1][j] > 0.0) 
       Value_1 += (Data[1][j] - Theory[1][j]) * (Data[1][j] - Theory[1][j]);

   Value_2 = 0.0;
   for( j=0; j<No_of_POINTS; j++ ) 
     if (Data[2][j] > 0.0) 
       Value_2 += (Data[2][j] - Theory[2][j]) * (Data[2][j] - Theory[2][j]);
   
   Value = 1.0/(double)No_of_POINTS * Value_0 + Value_1 + Value_2; 
   
   Value = sqrt(Value);

  }

  else Value = DBL_MAX; // 10.0; // Value = DBL_MAX;

  return(Value);
}

