#include <MODEL.h>

void Fixed_Points_nS1Cell(Parameter_Table * Table, 
                          double * Vector_Stationarity_Lower,
		                      double * Vector_Stationarity_Inter,
		                      double * Vector_Stationarity_Upper, 
                          double Epsilon );


void Fixed_Points_1S1Cell(Parameter_Table * Table, 
                          double * Vector_Stationarity_Lower,
		                      double * Vector_Stationarity_Inter,
		                      double * Vector_Stationarity_Upper );

void Fixed_Points_All( Parameter_Table * Table,
		       double * Vector_Stationarity_Lower,
		       double * Vector_Stationarity_Inter,
		       double * Vector_Stationarity_Upper,
		       double Epsilon)
{

  if(Table->No_of_RESOURCES == 1 && Table->No_of_CELLS == 1 )  
    Fixed_Points_1S1Cell(Table, Vector_Stationarity_Inter, Vector_Stationarity_Lower, Vector_Stationarity_Upper);
  else if (Table->No_of_CELLS == 1)
    Fixed_Points_nS1Cell(Table, Vector_Stationarity_Inter, Vector_Stationarity_Lower, Vector_Stationarity_Upper, Epsilon);
  else{
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics occurs in one single patch or local population\n");
    printf("However, No of CELLS is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  } 

}

void Fixed_Points_nS1Cell(Parameter_Table * Table, 
                          double * Vector_Stationarity_Lower,
		                      double * Vector_Stationarity_Inter,
		                      double * Vector_Stationarity_Upper, 
                          double Epsilon )
  
{ 
  /* This solution works for n species or types (with mutation) on one single local population */
  int i, k;
  double x, R_0; 
  double d, e, p, B, d_A;

  if(Table->No_of_CELLS > 1) {
      printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
      printf("Thefore, the dynamics occurs in one single patch or local population\n");
      printf("However, No of CELLS is larger than 1 (N = %d)\n", Table->No_of_CELLS);
      printf("The program will safely exit\n");
      Print_Press_Key(1,0,".");
      exit(0);
  }

  assert( Table->p_1 > 0.0 ); /* Mutation between strategies */

  double * Y = Vector_Stationary_Lower;
  double * y = (double *)calloc(Table->No_of_RESOURCES, sizeof(double)); 

  /* Define the tridiagonal system: */


  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  Vector_Stationarity_Lower[i];
      Vector_Stationarity_Upper[i] =                  Vector_Stationarity_Lower[i];
      Table->Vector_Model_Variables_Stationarity[i] = Vector_Stationarity_Lower[i];
  }

  free(y); 
}

void Fixed_Points_1S1Cell(Parameter_Table * Table, 
                          double * Vector_Stationarity_Lower,
		                      double * Vector_Stationarity_Inter,
		                      double * Vector_Stationarity_Upper )
  
{ /* This solution works for only one species or type (no mutation) on one single local population */
    int i, k;
    double x, R_0; 
    double d, e, p, B, d_A;

    assert( Table->p_1             == 0.0); /* No mutation because there is only one strategy */  
    assert_right_model_definition( Table );

    if(Table->No_of_RESOURCES > 0 || Table->No_of_CELLS > 1) {
      
   }

   p = Table->p_1; 
  
    for(i=0; i<Table->No_of_RESOURCES; i++) {
      B   = Table->Beta_AP[i]; 
      d_A = Table->Delta_AP[i];
      d   = Table->Delta_RP[i]; 
      e   = Table->Eta_RP[i]; 
      R_0 =  (1.0 - 2.0*p) * (B/d_A) / (1.0 + d/e);

      if (R_0 > 1.0) {
        k = i * 2 + 1; /* Adult Plant Density  m/N */
        Vector_Stationarity_Lower[k] = (R_0-1.0) * (1.0+d/e) / (R_0*(1.0+d/e) - 1);
        x = Vector_Stationarity_Lower[k];

        k = i * 2;      /* Propagule Density    n/N */
        Vector_Stationarity_Lower[k] =  x / (1.0 + x) * d_A/e;
      }
      else {
        k = i * 2 + 1; /* Adult Plant Density  m/N */
        Vector_Stationarity_Lower[k] = 0.0; 

        k = i * 2;      /* Propagule Density    n/N */
        Vector_Stationarity_Lower[k] =  0.0; 
      }   
    }
  
    for(i=0; i<Table->No_of_RESOURCES; i++) {
      k = i * 2 + 1; /* Adult Plant Number  m */
      Vector_Stationarity_Lower[k] *= Table->K_R; 

      k = i * 2;      /* Propagule Number n */
      x = Vector_Stationarity_Lower[k];
      Vector_Stationarity_Lower[k] *= Table->K_R;
    }
   
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  Vector_Stationarity_Lower[i];
      Vector_Stationarity_Upper[i] =                  Vector_Stationarity_Lower[i];
      Table->Vector_Model_Variables_Stationarity[i] = Vector_Stationarity_Lower[i];
    }
}

/* The following three functions are still under construction */
/* 
int  Coexistence_Condition ( Parameter_Table * Table )
{
  int Condition_Bool;
  
  return(Condition_Bool);
}

double Coexistence_Condition_Double ( Parameter_Table * Table )
{
  int Condition_Bool;
  double Condition_Double;

  
  Condition_Double = (double)Condition_Bool;
  return(Condition_Double);
}

double Function_to_Type_of_Stability( Parameter_Table * Table )
{
 
  *  Input arguments:
  *
  *  . Table, a pointer to the main data structure controling all parameters
  *  of the excution.

  *  Output arguments:
  *
  *  . Type_of_Stability, the value that takes this function for parameters values
  *  as defined in input Table.
 
  int i, k;
  int Type_of_Stability; 
  double Type_of_Stability_Double;

  
  Type_of_Stability_Double = (double)Type_of_Stability;
  return(Type_of_Stability_Double);
}
*/
