#include <MODEL.h>

extern gsl_rng * r; 

void Random_Error_Model_within_Boundaries_Model (Parameter_Model * Model,
						   Parameter_Space * S )
{
  int i, key;
  double value, lo_P, hi_P;

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);

    value = lo_P + gsl_rng_uniform(r) * (hi_P - lo_P);
    
    Vector_Entry_into_Error_Model_Model( value, key, Model);
  }
}

void Random_Error_Model_within_Boundaries_Table (Parameter_Table * Table,
						   Parameter_Space * S )
{
  int i, key;
  double value, lo_P, hi_P;

  for ( i=0; i < S->No_of_PARAMETERS; i++ ) {
    key = S->Parameter_Index[i];
    lo_P = gsl_vector_get(S->P_min, i);
    hi_P = gsl_vector_get(S->P_MAX, i);

    value = lo_P + gsl_rng_uniform(r) * (hi_P - lo_P);
    
    Vector_Entry_into_Error_Model_Table(value, key, Table);
  }
}




