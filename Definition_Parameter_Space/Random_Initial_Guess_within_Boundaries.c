#include <MODEL.h>

extern gsl_rng * r; 

void Random_Initial_Guess_within_Boundaries_Model (Parameter_Model * Initial_Guess,
						   Parameter_Space * S )
{
  int i, key;
  double value, lo_P, hi_P;

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);

    value = lo_P + gsl_rng_uniform(r) * (hi_P - lo_P);
    
    Vector_Entry_into_Parameter_Model(value, key, Initial_Guess );
  }
}

void Random_Initial_Guess_within_Boundaries_Table (Parameter_Table * Initial_Guess,
						   Parameter_Space * S )
{
  int i, key;
  double value, lo_P, hi_P;

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);

    value = lo_P + gsl_rng_uniform(r) * (hi_P - lo_P);
    
    AssignVectorEntry_to_Structure(Initial_Guess, key, value);
  }
}




