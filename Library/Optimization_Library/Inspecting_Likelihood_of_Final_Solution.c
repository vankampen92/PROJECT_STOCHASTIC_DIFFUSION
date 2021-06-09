#include <MODEL.h>

int Model_Two_Phases( Parameter_Table * , Parameter_Fitting * ); 

double Inspecting_Likelihood_of_Final_Solution( const gsl_vector * x, void * Par )
{
  /* This GSL function allows an evaluation of function to minimize

     Output value:
     . Value
  */
  double Value, Value_j;
  int i,key;

  Parameter_Fitting * F = (Parameter_Fitting *)Par;

  int No_of_POINTS         = F->Data->No_of_POINTS;
  int No_of_VARIABLES      = F->Data->No_of_VARIABLES;
  double ** Data           = F->Data->N;
  double ** Theory         = F->Table->Matrix_Output_Variables;
  Parameter_Table * Table  = F->Table;
  Parameter_Space * Space  = F->Space;
  int No_of_PARAMETERS     = F->Space->No_of_PARAMETERS;
  int x_is_BOUNDED;
  int P_0, P_1, P_2;

  x_is_BOUNDED = Checking_for_Parameter_Boundaries( F, x );

  P_0 = Space->No_of_PARAMETERS;
  P_1 = Table->No_of_IC;
  P_2 = Table->No_of_ERROR_PARAMETERS;

  if (F->Verbose == 1) {

    if (x_is_BOUNDED == 1) {
      printf("The final solution is bounded:\n");

      for( i = 0;  i < P_0; i++ ) {
	key = Space->Parameter_Index[i];
	printf("%s=%10.3e  ", Table->Symbol_Parameters[key], gsl_vector_get (x, i) );
      }
      for( i = P_0; i < (P_0+P_1) ; i++ ) {
	key = Table->IC_Space->Parameter_Index[i-P_0];
	printf("%s=%10.3e  ", Table->Model_Variable_Symbol[key], gsl_vector_get (x, i) );
      }
      for( i = (P_0+P_1); i < (P_0+P_1+P_2); i++ ) {
	key = Table->E_Space->Parameter_Index[i-P_0-P_1];
	printf("Error =%10.3e  ", gsl_vector_get (x, i) );
      }
      printf("\n");
    }
    else printf("The final solution is NOT bounded (so, not written) :\n");

  }

  if(No_of_PARAMETERS > 0)
      Vector_Entries_into_Parameter_Table ( x, Table,
					    Space->Parameter_Index, No_of_PARAMETERS );

  if(Table->No_of_IC > 0)
      Vector_Entries_into_Parameter_Table_Initial_Condition ( x, Table,
							      Table->IC_Space->Parameter_Index,
							      No_of_PARAMETERS,
							      Table->No_of_IC );
  if(Table->No_of_ERROR_PARAMETERS > 0)
      Vector_Entries_into_Parameter_Table_Error_Model ( x, Table,
							Table->E_Space->Parameter_Index,
							No_of_PARAMETERS,
							Table->No_of_IC,
							Table->No_of_ERROR_PARAMETERS );
  int State;

  /* These two lines should be activited in case we require a    */
  /* 2-phase numerical integration (from 2000 to 2016 and from   */
  /* 2016 to 2033) of the full ODE system. This is sometimes     */
  /* needed if projected trajectories though extrapolation are   */
  /* to be calculated.                                           */
  /* if (F->TWO_PHASES == 1) State = Model_Two_Phases(Table, F); */
  /* else                    State = M_O_D_E_L(Table);           */

  State = M_O_D_E_L(Table); 

  assert(State == GSL_SUCCESS);

  Value = 0.0;
  for( i=0; i<No_of_VARIABLES; i++ ) {
    key = Table->OUTPUT_VARIABLE_INDEX[i];

    if (F->Verbose == 1)
      printf(" %s (%s):\n",
	     Table->Output_Variable_Name[key], Table->Output_Variable_Symbol[key]);
    
    Value_j = GSL_neglog_Error_Probability_Model( Data[i], Theory[i],
						  No_of_POINTS, No_of_VARIABLES, F,
						  GSL_neglog_Error_Probability_Model_Gaussian );
    Value += Value_j;
  }

  if( x_is_BOUNDED == 1 ) {
    if (F->Verbose == 1) {
      printf("Under Gaussian distributed independent errors, NLL (Data|Theory) = %g\n",
	     Value);
       #if defined CPGPLOT_REPRESENTATION
       GSL_CPGPLOT_Minimization_Simplex (F, x, 0, F->Function);
       #endif
    }
  }
  else {
    if (F->Verbose == 1)  {
      printf(" Since the final solution is NOT bounded!!!\n");
      printf(" Arbitrarily, we set the likelihood to the upper maximum value:\n");
      printf(" NLL (Data|Theory) = %g\n", Value);
    }

    Value = DBL_MAX;
  }

  return(Value);
}
