#include <MODEL.h>

double definition_Scan_Function( Parameter_Table * P )
{
  double x;
  double * y;

  int j = P->OUTPUT_VARIABLE_INDEX[0];

  x = definition_OutPut_Variables( j, y, 0.0, P );

  return (x);
}

double definition_OutPut_Variables(int j, double * y, const double t, Parameter_Table * Table)
{
  int i;
  double x = 0.0;

  /* Definition of the numerical order, from 0 to n, of model variables */
  #include <Model_Variables_Code.Include.c>

  /* Genuine output variables are those that derive from the state of the system and,
   * therefore, should be evaluated as a funcion of system state variables
   */

  if (j < Table->No_of_RESOURCES ) {
    Table->Focal_Resource = j;
    x = Total_Population_Resource_Species (y, Table);
  }
  else if (j < Table->OUTPUT_VARIABLES_GENUINE) {
    /* Derived output variables from model dynamic variables and parameters */
    j -= Table->No_of_RESOURCES;
   
    switch(j)
    {
    
    case 0:  x = Average_Individual_Density ( y, Table );
      break;
    
    case 1:  x = Standard_Deviation_Density ( y, Table );
      break;

    case 2:  x = Total_Population(y, Table);
      break;
      
    default:
      printf(" INVALID VARIABLE CODE in ./definition_OutPut_Variables.c (variable = %d)\n", j);
      printf(" Output Variables range from 0 to 2: Program will exit!  ");
      printf(" Press any key..."); getchar();
      exit(0);
    }
  }
  /* The last MODEL_STATE_VARIABLES output variables are the MODEL_STATE_VARIABLES */
  else {
    j -= Table->OUTPUT_VARIABLES_GENUINE; /* #defined in MODEL.h */
    assert( j < Table->MODEL_STATE_VARIABLES );
    x = y[j];
  }

  return (x);
}
