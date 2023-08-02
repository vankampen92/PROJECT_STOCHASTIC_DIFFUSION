#include <MODEL.h>

double neglogBinomial_Free_Consumers (double , Parameter_Table * );

double GSL_Function_to_Minimize_Binomial_Free_Consumers( const gsl_vector * x, void * Par )
{
  /* This GSL function allows the evaluation of a likelihood function to minimize

     This likelihood corresponds to the probability of observing a number of free 
     consumers at stationarity for the Holling Type II model that isolates 
     feeding dynamics. This is the reason there is an assert to verify the 
     corret model has been chosen and compiled. 

     Output value:
     . Value
  */
  double Value;
  int i,j;

  Parameter_Fitting * F   = (Parameter_Fitting *)Par;

  int No_of_POINTS         = F->Data->No_of_POINTS;
  int No_of_VARIABLES      = F->Data->No_of_VARIABLES;
  double ** Data           = F->Data->N;
  double ** Theory;
  Parameter_Table * Table  = F->Table; 
  Parameter_Space * Space  = F->Space;
  int No_of_PARAMETERS     = F->Space->No_of_PARAMETERS;
  int x_is_BOUNDED; 
                                                                      /* 12: DIFFUSION_HII_1D */
  assert ( Table->TYPE_of_MODEL == 9 || Table->TYPE_of_MODEL == 12);  /*  9: DIFFUSION_HII_2D */
  assert ( Table->SUB_OUTPUT_VARIABLES == 1 );                        /* Only one!!!      */
  assert ( No_of_VARIABLES == 1 );                                    /* Only one!!!      */

  Theory = (double **)calloc(Table->SUB_OUTPUT_VARIABLES, sizeof(double *) ); 
  for(i=0; i<Table->SUB_OUTPUT_VARIABLES; i++)
    Theory[i] = (double *)calloc( No_of_POINTS, sizeof(double) ); 
    
  x_is_BOUNDED = 1; /* By default, the neg Log Likelihood is evaluated */
  
  if(F->Minimization == 1)
    x_is_BOUNDED = Checking_for_Parameter_Boundaries( F, x );
  
  if( x_is_BOUNDED == 1 ) { 
						    
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
    
    for(i=0; i<Table->SUB_OUTPUT_VARIABLES; i++)
      for(j=0; j<No_of_POINTS; j++)
	      /* Probability of data given a choice of parameter values: 
          -log( P(Data | parameters) ) 
        */ 
	      Theory[i][j] = neglogBinomial_Free_Consumers (Data[i][j], Table);
      
    int Theory_is_NOT_a_NUMBER = 0;
    for( i=0; i<No_of_VARIABLES; i++ )
      Theory_is_NOT_a_NUMBER += da_vector_isnan(Theory[i], No_of_POINTS);

    if( Theory_is_NOT_a_NUMBER == 0 ) { 
      Value = 0.0;
      for( i=0; i<No_of_VARIABLES; i++ )
	      for(j=0; j<No_of_POINTS; j++)
	        Value += Theory[i][j]; 
    }
    else Value = DBL_MAX; /* Likelihood takes a non-feasible value */
    
  }
  
  else Value = DBL_MAX;   /* Parameter combination outside boundaries */

  for( i=0; i<No_of_VARIABLES; i++ ) free(Theory[i]); 
  free(Theory);
  
  return(Value);
}

double neglogBinomial_Free_Consumers (double Data, Parameter_Table * Table)
{
  double K_R, y_R, Theta;
  double Theory;
  double p;
  int    k, n; 

  K_R   = (double)Table->K_R;
  y_R   = Table->TOTAL_No_of_RESOURCES; /* Resource density (in number of resource units) */
  Theta = Table->Alpha_C_0 * y_R/K_R;
  
  p = Table->Nu_C_0 / (Table->Nu_C_0 + Theta) ;
  k = (unsigned int)Data;
  n = Table->TOTAL_No_of_CONSUMERS;

  Theory = gsl_ran_binomial_pdf(k, p, n); 

  Theory = -log(Theory); 
  
  return(Theory);
}
