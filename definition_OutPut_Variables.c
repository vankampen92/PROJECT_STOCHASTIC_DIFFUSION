#include <MODEL.h>

double definition_OutPut_Variables(int j, double * y, const double t, Parameter_Table * Table)
{
  int i, N;
  double x = 0.0;

  /* Definition of the numerical order, from 0 to n, of model variables */
  #include <Model_Variables_Code.Include.c>

  /* Genuine output variables are those that derive from the state of the system and,
   * therefore, should be evaluated as a funcion of system state variables
   */

  if (Table->T->TYPE_of_TIME_DEPENDENCE == 0 ) N = 0; 
  else                                         N = Table->TDC->TIME_DEPENDENT_PARAMETERS;

  if (j < Table->LOCAL_STATE_VARIABLES) {
    Table->Focal_Resource = j;
    x = Total_Population_Resource_Species (y, Table);
    /* These are Global State Variables, this is, they are simply the sum of 
       Total Populations across the entire Spatial Metapopulation for every 
       species type. In other words, GLOBAL STATE VARIABLES are the LOCAL 
       STATE VARIABLES summed over the whole network of interconnected local 
       communities. 
       
       Of course, the number of global variables is equal to then number of 
       LOCAL_STATE_VARIABLES. 
    */
  }
  else if (j < (Table->OUTPUT_VARIABLES_GENUINE) ) {
    /* Recall what Table->OUTPUT_VARIABLES_GENUINE stores: 

        OUTPUT_VARIABLES_GENUINE = LOCAL_STATE_VARIABLES + OUTPUT_VARIABLES_TRUE_DERIVED
    
    */
    /* Derived output variables from model dynamic variables and parameters */
    j -= Table->LOCAL_STATE_VARIABLES;
   
    switch(j)
    {
    
    case 0:  x = Average_Individual_Density ( y, Table );
      break;
    
    case 1:  x = Standard_Deviation_Density ( y, Table );
      break;

    case 2:  x = Total_Population(y, Table);
      break;

    case 3:  x = Total_Population_Resources(y, Table);
      break;
    
    case 4:  x = Total_Population_Free_Consumers ( y, Table );
      break;

    case 5:  x = Total_Population_Mature_Consumers(y, Table);
      break;

    case 6:  x = Total_Population_Handling_Consumers ( y, Table );
      break;

    case 7:  x = Total_Population_Triplet_Consumers(y, Table);
      break;

    case 8:  x = Total_Population_Consumers(y, Table);
      break;
    
    default:
      printf(" INVALID VARIABLE CODE in ./definition_OutPut_Variables.c (variable = %d)\n", j);
      printf(" Output Variables range from 0 to 2: Program will exit!  ");
      printf(" Press any key..."); getchar();
      exit(0);
    }
  }
  /* The last output variables are the MODEL_STATE_VARIABLES */
  else {
    j -= Table->OUTPUT_VARIABLES_GENUINE; /* #defined in MODEL.h */
    assert( j < Table->MODEL_STATE_VARIABLES );
    x = y[j];
  }

  return (x);
}

double definition_Scan_Function( Parameter_Table * P )
{
  double x;
  double * y;

  int j = P->OUTPUT_VARIABLE_INDEX[0];

  x = definition_OutPut_Variables( j, y, 0.0, P );

  return (x);
}
