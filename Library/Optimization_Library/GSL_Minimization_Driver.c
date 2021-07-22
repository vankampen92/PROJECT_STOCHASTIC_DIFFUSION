#include <MODEL.h>

extern int MY_ERROR_HANDLER;

double GSL_Minimization_Driver( Parameter_Fitting * F )
{
  /* This function is just a wrapper for a GSL minimization algorithm:

        GSL_Minimization_Simplex( Parameter_Fitting * F,
				  gsl_vector * Initial_Guess,
				  gsl_vector * Solution,
				  double ( * Function) (const gsl_vector *, void * ) )
     Input arguments:

     Output arguments:

  */
  int i, j, k;
  
  double Value; 
  
  Parameter_Table * Table = F->Table;
  Parameter_Space * Space = F->Space; 
  
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
  
  MY_ERROR_HANDLER = 0;
  gsl_error_handler_t * old_handler = gsl_set_error_handler ( &my_error_handler );

  if( F->Minimization == 1 ) Value = GSL_Minimization_Simplex( F, x, x, F->Function );
  
  else if ( F->Minimization == 0 ) {
    Value =  ( * F->Function )( x, F );
    #if defined VERBOSE
    #if defined CPGPLOT_REPRESENTATION
       GSL_CPGPLOT_Minimization_Simplex (F, x, 0, F->Function);
    #endif
    #endif 
    return(Value);
  }
  
  else     printf(" Error in 1/0 Minimization input argument!\n ---> Minimization = %d\n",
           F->Minimization );

  gsl_set_error_handler (old_handler);

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

  if ( F->Minimization == 0 ) 
    Value = Inspecting_Likelihood_of_Final_Solution(x, F);

  if (MY_ERROR_HANDLER == 1) {
    Value = DBL_MAX;
    MY_ERROR_HANDLER = 0;
  }

  F->Bounded_Parameter_Set = Checking_for_Parameter_Boundaries(F, x);
  
  gsl_vector_free(x);
  
  return(Value); 
}
