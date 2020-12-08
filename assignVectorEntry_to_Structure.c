#include <MODEL.h>

void Vector_Entries_into_Parameter_Table ( const gsl_vector * X, Parameter_Table * P,
					   int * Parameter_Index, int No_of_PARAMETERS )
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_PARAMETERS; i++) {
    key = Parameter_Index[i];
    value = gsl_vector_get(X, i);
    AssignVectorEntry_to_Structure( P, key, value );
  }
}

void Vector_Entries_into_Parameter_Table_Initial_Condition ( const gsl_vector * X,
							     Parameter_Table * P,
							     int * Index,
							     int No_of_PARAMETERS,
							     int No_of_IC)
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_IC; i++) {
    key = Index[i];
    value = gsl_vector_get(X, i+No_of_PARAMETERS);
    Vector_Entry_into_Model_Variable_Initial_Condition_Table( value, key, P );
  }
}

void AssignVectorEntry_to_Structure(Parameter_Table * P, int j, double value)
{

  switch(j)
    {
   
    case  0: P->Mu = value;     
      break;
    case  1: P->No_of_INDIVIDUALS = value;  
      break;
    case  2: P->No_of_CELLS = value;  
      break;
    case  3: P->No_of_CELLS_X = value;  
      break;
    case  4: P->No_of_CELLS_Y = value;  
      break;
    case  5: P->No_of_SPECIES;  
      break;
    
    default:
      printf(".... INVALID PARAMETER KEY (key = %d)\n", j);

      printf(".... The permited correspondences are:\n");
      printf("\n");
      fprintf_Model_Parameters(stdout, P);

      printf(" The maximum number of parameters is Number_PAR_MAX\n");
      printf(" The permited number of keys go from 0, to %d\n", MODEL_PARAMETERS_MAXIMUM-1);

      exit(0);
    }
}
