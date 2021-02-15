#include <MODEL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int function (double t, const double y[], double dydt[], void *params)
{
  int i, n, j;
  int Sp; 
  
  Parameter_Table * Table = (Parameter_Table *)params;

  Sp = Table->No_of_RESOURCES; 

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>

  Resetting_Lambda_Delta_Vectors (Table); 

  n= 0; 
  for (j=0; j<Table->No_of_CELLS; j++) 
    for(i=0; i<Sp; i++) {
      n = j*Sp + i; 
      dydt[n] = Table->Lambda_R[i] - Table->Delta_R[i] * y[n];
      n++;
    }
       
  n= 0; 
  for (j=0; j<Table->No_of_CELLS; j++) { 
   
    for(i=0; i<Sp; i++) { 
      dydt[n] += In_Mu(Table, n, i, j, y) - Out_Mu_Per_Capita(Table, i, j) * y[n];;
      n++;
    }
  }
  
  return GSL_SUCCESS;
}
