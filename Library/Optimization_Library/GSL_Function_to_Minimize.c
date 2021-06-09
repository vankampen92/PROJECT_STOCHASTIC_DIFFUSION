#include <MODEL.h>

double GSL_Function_to_Minimize( const gsl_vector * x, void * Par )
{
  /* This GSL function allows an evaluation of function to minimize

     Output value:
     . Value
  */
  double Value,t;
  int i,j,k;

  Parameter_Fitting * F = (Parameter_Fitting *)Par;

  int No_of_POINTS         = F->Data->No_of_POINTS;
  int No_of_VARIABLES      = F->Data->No_of_VARIABLES;
  double ** Data           = F->Data->N;
  double ** Theory         = F->Table->Matrix_Output_Variables;
  Parameter_Table * Table  = F->Table;
  Parameter_Space * Space  = F->Space;   
  int No_of_PARAMETERS     = F->Space->No_of_PARAMETERS;
  int x_is_BOUNDED;
  int Theory_is_NOT_a_NUMBER; 
  
  x_is_BOUNDED = Checking_for_Parameter_Boundaries( F, x );

  if( x_is_BOUNDED == 1 ) {

    Vector_Entries_into_Parameter_Table ( x, Table,
					  Space->Parameter_Index, No_of_PARAMETERS );

    if(Table->No_of_IC > 0)
      Vector_Entries_into_Parameter_Table_Initial_Condition ( x, Table,
							      Table->IC_Space->Parameter_Index,
							      No_of_PARAMETERS,
							      Table->IC_Space->No_of_PARAMETERS);

    if ( Table->T->TYPE_of_TIME_DEPENDENCE > 0 ) {

      int TYPE_0_PARAMETERS  = Table->TDC->TYPE_0_PARAMETERS;
      int TYPE_1_PARAMETERS  = Table->TDC->TYPE_1_PARAMETERS;
      int TYPE_2_PARAMETERS  = Table->TDC->TYPE_2_PARAMETERS;

      if ( !Table->x_Bool ) {
        for(i = 0; i < TYPE_2_PARAMETERS; i++) {
          k = i+TYPE_0_PARAMETERS+TYPE_1_PARAMETERS;
          for(j = 0; j<Table->TDC->No_of_TIMES; j++) {
            t = Table->T->Time_Vector[j];
            Table->TDC->Dependent_Parameter[k][j]=Time_Dependence_Resolve(Table, Table->TDC->Index_Dependent_Parameters[k], Table->TDC->Forcing_Pattern_Parameters[k], t);
          }
        }
      }
    }

    int State = M_O_D_E_L(Table);

    Theory_is_NOT_a_NUMBER = 0;
    for( i=0; i<No_of_VARIABLES; i++ )
      Theory_is_NOT_a_NUMBER += da_vector_isnan(Theory[i], No_of_POINTS);
  
    if( Theory_is_NOT_a_NUMBER == 0 && State == GSL_SUCCESS) { 
      Value = 0.0;
      for( i=0; i<No_of_VARIABLES; i++ )
	for( j=0; j<No_of_POINTS; j++ )
	  Value += (Data[i][j] - Theory[i][j]) * (Data[i][j] - Theory[i][j]) / Theory[i][j];

      Value = sqrt(Value);
    }
    else Value = DBL_MAX;
    
  }
  else Value = DBL_MAX;

  return(Value);
}
