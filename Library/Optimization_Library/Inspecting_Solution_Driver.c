#include <MODEL.h>

double Inspecting_Solution_Driver( Parameter_Fitting * F )
{
  /* This function is just a wrapper for:

     double Inspecting_Likelihood_of_Final_Solution (Parameter_Fitting *)

     Input arguments:
     . Parameter_Fitting * F

     Output arguments:
     . Value, likelihood of model parameters passed on the function
       through the F Parameter_Fitting structure.
     . Graphical representation of the solution

     Both outputs are obtained through calling:

     Inspecting_Likelihood_of_Final_Solution( gsl_vector * , Parameter_Fitting *)

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

  Value = Inspecting_Likelihood_of_Final_Solution(x, F);

  F->Bounded_Parameter_Set = Checking_for_Parameter_Boundaries(F, x);

  gsl_vector_free(x);

  return(Value);
}
