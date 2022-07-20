#include <MODEL.h>

double neglogStarionary_Ditribution_Beddington_DeAngelis(int , int,
							 Parameter_Table * );

double GSL_Function_to_Minimize_Beddington_DeAngelis( const gsl_vector * x, void * Par )
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
  int n, m;

  Parameter_Fitting * F   = (Parameter_Fitting *)Par;

  int No_of_POINTS         = F->Data->No_of_POINTS;
  int No_of_VARIABLES      = F->Data->No_of_VARIABLES;
  double ** Data           = F->Data->N;
  double * Theory;
  Parameter_Table * Table  = F->Table; 
  Parameter_Space * Space  = F->Space;
  int No_of_PARAMETERS     = F->Space->No_of_PARAMETERS;
  int x_is_BOUNDED; 
                                                                      
  assert ( Table->TYPE_of_MODEL == 13 );                              /* 13: DIFFUSION_BD_2D */
  
  Theory = (double *)calloc( No_of_POINTS, sizeof(double) ); 
    
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
    // for(i=0; i < No_of_VARIABLES; i++) {
    for(j=0; j<No_of_POINTS; j++) {
      
	n = (int)Data[0][j];
	m = (int)Data[1][j];
    
	/* Probability of each data across realization: -log( P(Data | parameters) ) */ 
	Theory[j] = neglogStarionary_Ditribution_Beddington_DeAngelis (n, m, Table);
    }
    // }
  
    int Theory_is_NOT_a_NUMBER = da_vector_isnan(Theory, No_of_POINTS);

    if( Theory_is_NOT_a_NUMBER == 0 ) { 
      Value = 0.0;
      
      for(j=0; j<No_of_POINTS; j++)
	  Value += Theory[j]; 
    }
    else Value = DBL_MAX; /* Likelihood takes a non-feasible value */   

  }
  
  else Value = DBL_MAX;   /* Parameter combination outside boundaries */
  
  free(Theory);
  
  return(Value);
}

double neglogStarionary_Ditribution_Beddington_DeAngelis(int n,
							 int m,
							 Parameter_Table * Table)
{
  double Theory;
  
  Stationary_Probability_Distribution (Table);                /* True Calculation    */

  assert(Table->MEq->n_DIMENSION == 2); 
  
  Theory = -log( Table->MEq->PS_nm[n][m] );
  
  return(Theory);
}
