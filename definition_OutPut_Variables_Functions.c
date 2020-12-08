#include <MODEL.h>

double Total_Population(double * Y, Parameter_Table * Table)
{
  double x;
  int i;
  
  if (Table->TYPE_of_MODEL == 0) {
    x = 0.0;
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) x += Y[i]; 
  }

  return(x); 
}


double Average_Individual_Density  ( double * y, Parameter_Table * Table )
{
  double x, x_0, x_1;
  int i;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
 
   if (Table->TYPE_of_MODEL == 0) {
      x = Average_double_Vector(y, Table->No_of_CELLS);
   }
   else {
     printf(" Type of Model is ill-defined. Check input arguments. TYPE of MODEL is %d\n",
	    Table->TYPE_of_MODEL);
     printf(" The program will exit\nm");
     exit(0);
   }
   
  return (x);
}

double Standard_Deviation_Density  ( double * y, Parameter_Table * Table )
{
  double x, xm;
  int i;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
 
   if (Table->TYPE_of_MODEL == 0) {
     xm = Average_double_Vector(y, Table->No_of_CELLS);
     x  = Variance_double_Vector(y, Table->No_of_CELLS);
     x = sqrt( x - xm*xm );  
   }
   else {
     printf(" Type of Model is ill-defined. Check input arguments. TYPE of MODEL is %d\n",
	    Table->TYPE_of_MODEL);
     printf(" The program will exit\nm");
     exit(0);
   }
   
  return (x);
}

