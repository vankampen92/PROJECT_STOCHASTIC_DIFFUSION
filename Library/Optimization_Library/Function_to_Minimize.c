#include <MODEL.h>

double Function_to_Minimize( Parameter_Table * Table )
{
  /* This function is just a wrapper for a generic 
   *				  
   *	double ( * Function) (const gsl_vector *, void * ) )

   * All the 'parafernalia' is required since parameter values should be
   * conveniently written onto a gsl_vector x as prerequisite to 
   * properly call the function 
     	
   *  Input arguments:
   *
   *  . Table, a pointer to the main data structure controling all parameters
   *  of the excution.
  
   *  Output arguments:
   *
   *  . value, the value that takes this function for parameters values 
   *  as defined in input Table. 
   */
  
  int i, j, k;
  
  double Value; 
 
  Parameter_Fitting * F   = (Parameter_Fitting *)Table->Fitting_Data;
  Parameter_Space * Space  = F->Space;   
  
  gsl_vector * x  = gsl_vector_alloc(Table->S->No_of_PARAMETERS
				     + Table->No_of_IC + Table->No_of_ERROR_PARAMETERS);

  if(Space->No_of_PARAMETERS > 0) 
    Parameter_Table_into_Vector_Entries ( Table, x,
					  Table->S->Parameter_Index,
					  Table->S->No_of_PARAMETERS );

  assert(Space->No_of_PARAMETERS > 0); 
  
  if(Table->No_of_IC > 0) 
    Parameter_Table_into_Vector_Entries_Initial_Condition ( Table, x,
							    Table->IC_Space->Parameter_Index,
							    Table->S->No_of_PARAMETERS,
							    Table->IC_Space->No_of_PARAMETERS);
  if(Table->No_of_ERROR_PARAMETERS > 0)
    Parameter_Table_into_Vector_Entries_Error_Model ( Table, x,
						      Table->E_Space->Parameter_Index,
						      Table->S->No_of_PARAMETERS,
						      Table->No_of_IC,
						      Table->No_of_ERROR_PARAMETERS);
  
  Value =  ( * F->Function )( x, F );
   
  if(Space->No_of_PARAMETERS > 0) 
    Vector_Entries_into_Parameter_Table ( x, Table,
					Table->S->Parameter_Index, Table->S->No_of_PARAMETERS );

  if(Table->No_of_IC > 0)
    Vector_Entries_into_Parameter_Table_Initial_Condition ( x, Table, 
							    Table->IC_Space->Parameter_Index,
							    Table->S->No_of_PARAMETERS,
							    Table->IC_Space->No_of_PARAMETERS);
  if(Table->No_of_ERROR_PARAMETERS > 0)
    Vector_Entries_into_Parameter_Table_Error_Model ( x, Table, 
						      Table->E_Space->Parameter_Index,
						      Table->S->No_of_PARAMETERS,
						      Table->No_of_IC,
						      Table->No_of_ERROR_PARAMETERS );

  if ( Table->T->TYPE_of_TIME_DEPENDENCE > 0 ) {

    int TYPE_0_PARAMETERS  = Table->TDC->TYPE_0_PARAMETERS;
    int TYPE_1_PARAMETERS  = Table->TDC->TYPE_1_PARAMETERS;
    int TYPE_2_PARAMETERS  = Table->TDC->TYPE_2_PARAMETERS;

    if ( !Table->x_Bool ) {
      for(i = 0; i < TYPE_2_PARAMETERS; i++) {
        k = i+TYPE_0_PARAMETERS+TYPE_1_PARAMETERS;
        for(j = 0; j<Table->TDC->No_of_TIMES; j++) {
          Table->TDC->Dependent_Parameter[k][j]=Time_Dependence_Resolve(Table, Table->TDC->Index_Dependent_Parameters[k], Table->TDC->Forcing_Pattern_Parameters[k], Table->T->Time_Vector[j]);
        }
      }
    }
  }

  gsl_vector_free(x);
  
  return(Value); 
}
