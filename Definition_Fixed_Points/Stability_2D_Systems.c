#include <MODEL.h>
                            
double Calculate_Type_of_Stability_2D( Parameter_Table * Table)
{
 
  /*  Input arguments:
  *
  *  . Table, a pointer to the main data structure controling all model 
  *    and execution parameters

  *  Output arguments:
  *
  *  . Type_of_Stability, the value that takes this function for parameters values
  *    as defined in input Table. It is a double number that can be used to
  *    plot the type of stability of the system across the parameter space.
  */
  int Type_of_Stability;
  double Tra, Det, Dis;
  double Type_of_Stability_Double;

  double * Y = Table->Vector_Model_Variables_Stationarity;

  assert(Table->MODEL_STATE_VARIABLES == 2);

  /* BEGIN: Stability Analysis */ 
  Tra = Trace(Table, Y); 
  Det = Determinant(Table, Y);
  Dis = Tra * Tra - 4.0 * Det;
   
  if (Det < 0.0) {
    /* The system is unstable */
    Type_of_Stability = 0; /* There is always one positive and one negative eigen value */
  }
  else {
    /* */
    if (Tra < 0.0) {
      if (Dis > 0.0) {
        Type_of_Stability = 2;  /* Stable Node */
      }
      else {
        Type_of_Stability = 3;  /* Stable Focus */
      }
    }
    else {
      if (Dis > 0.0) {
        Type_of_Stability = 0;  /* Unstable Node */
      }
      else {
        Type_of_Stability = 1;  /* Unstable Focus */
      }
    }
  }
  /*   END: Stability Analysis */ 
  
  Type_of_Stability_Double = (double)Type_of_Stability;

  return(Type_of_Stability_Double);
}

double Calculate_Turing_Instabilities( Parameter_Table * Table)
{
 
  /*  Input arguments:
  *
  *  . Table, a pointer to the main data structure controling all model 
  *    and execution parameters

  *  Output arguments:
  *
  *  . Turing_Structures, the booleanm value that takes this function for the 
  *    parameters values as defined in input Table. It is a double number that 
  *    it 0 if there are not Turing instabilities and 1 if they are. It can be 
  *    used to plot the exsitence of Turing instabilities in the reaction-diffusion 
  *    system across the parameter space.
  */

  double Tra, Det, Tra_D, d;
  double Turing_Structures;

  double * Y = Table->Vector_Model_Variables_Stationarity;

  assert(Table->MODEL_STATE_VARIABLES == 2);

  /* BEGIN: Stability Analysis */ 
  Tra = Trace(Table, Y); 
  Det = Determinant(Table, Y);

  d = Table->Mu_RP[1]/Table->Mu_RP[0]; 
   
  if (Det > 0.0 && Tra < 0.0) {
    /* The system may have Turing Instabilities */
    Tra_D = Trace_Diffusion_Matrix(Table, Y);

    if( Tra_D > (2.0  * sqrt(Det) * sqrt(d)) ) /* Turing Condition */  
    {
      /* The system does have Turing Instabilities */
      Turing_Structures = 1.0;
    }
    else {
      Turing_Structures = 0.0;
    }
  }
  else {
    Turing_Structures = 0.0;
  }

  return(Turing_Structures);
}

double Trace_Diffusion_Matrix( Parameter_Table * Table, double * Y)
{
  double Tra, d;

  assert(Table->MODEL_STATE_VARIABLES == 2);
  gsl_matrix * m = gsl_matrix_alloc( Table->MODEL_STATE_VARIABLES, 
                                     Table->MODEL_STATE_VARIABLES); /* allocating matrix m */
  
  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  /* Optimally, this Function should be inline... */
  JACOBIAN_Matrix(m, Y, 0.0, Table->K, Table);
  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  
  d = Table->Mu_RP[1]/Table->Mu_RP[0]; 

  Tra = d * gsl_matrix_get(m, 0, 0) + gsl_matrix_get(m, 1, 1);

  gsl_matrix_free(m); /* releasing memory */
  return(Tra);
}

double Trace( Parameter_Table * Table, double * Y)
{
  double Tra; 

  assert(Table->MODEL_STATE_VARIABLES == 2);
  gsl_matrix * m = gsl_matrix_alloc( Table->MODEL_STATE_VARIABLES, 
                                     Table->MODEL_STATE_VARIABLES); /* allocating matrix m */

  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  /* Optimally, this Function should be inline... */
  JACOBIAN_Matrix(m, Y, 0.0, Table->K, Table);
  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  
  Tra = gsl_matrix_get(m, 0, 0) + gsl_matrix_get(m, 1, 1);

  gsl_matrix_free(m); /* releasing memory */
  return(Tra);
}

double Determinant( Parameter_Table * Table, double * Y)
{
  double Det; 

  assert(Table->MODEL_STATE_VARIABLES == 2);
  gsl_matrix * m = gsl_matrix_alloc( Table->MODEL_STATE_VARIABLES, 
                                     Table->MODEL_STATE_VARIABLES); /* allocating matrix m */

  /* Setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  /* Optimally, this Function should be inline... */
  JACOBIAN_Matrix(m, Y, 0.0, Table->K, Table);
  /* End of setting the Jacobian matrix evaluated at (y[0], ..., y[K]) */
  
  Det = gsl_matrix_get(m, 0, 0) * gsl_matrix_get(m, 1, 1) - gsl_matrix_get(m, 0, 1) * gsl_matrix_get(m, 1, 0);

  gsl_matrix_free(m); /* releasing memory */
  return(Det);
}

