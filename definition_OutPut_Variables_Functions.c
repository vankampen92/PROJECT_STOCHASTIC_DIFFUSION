#include <MODEL.h>

double Local_Population_Resources(int i, const double * Y, Parameter_Table * Table)
{
  /* Input: 
  
      . i : the i-th Local Population or Patch
      . Y : the total state vector 
      . Table:  Full Table of Parameters 
    
     Output:

      . y_S: total local population of resources over all resources types or species
        at the i-th local patch. 
  */
  double y_S; 
  int k, n; 

  assert(Table->TYPE_of_MODEL == 20 || Table->TYPE_of_MODEL == 21); 

  y_S = 0.0;
  for(k = 0; k < Table->No_of_RESOURCES; k++) {
      if(Table->TYPE_of_MODEL == 20)
        n   = i*Table->LOCAL_STATE_VARIABLES + 2*k+1;
      else if (Table->TYPE_of_MODEL == 21)
        n   = i*Table->LOCAL_STATE_VARIABLES + k;
      else {
        printf(" Only Table->TYPE_of_MODEL == 20 and Table->TYPE_of_MODEL == 21  can make use\n"); 
        printf(" of the Local_Population_Resources (...); (definition_OutPut_Variables_Functions.c)\n");
        printf(" at this state\n");
        assert(Table->TYPE_of_MODEL == 20 || Table->TYPE_of_MODEL == 21); 
      }       

      y_S += Y[n];
  }

  return (y_S);
}

double Total_Population_Resources(double * Y, Parameter_Table * Table)
{
  double x;
  int i;

  /* Free Consumers
  */
  x = 0.0;

  assert(Table->TYPE_of_MODEL == 3); 
  
  if (Table->TYPE_of_MODEL == 3) {
    x = 0.0;
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
      if(i%Table->LOCAL_STATE_VARIABLES == 0) x += Y[i]; 
    
  }

  return(x); 
}

double Total_Population_Free_Consumers(double * Y, Parameter_Table * Table)
{
  double x;
  int i,j;

  /* Free Consumers
  */
  x = 0.0;

  assert(Table->TYPE_of_MODEL == 3); 
  
  if (Table->TYPE_of_MODEL == 3) {
    x = 0.0;
    for(j=0; j<Table->N_E; j++){
      for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
	      if(i%Table->LOCAL_STATE_VARIABLES == (j+1)) x += Y[i]; 
    }
  }

  return(x); 
}

double Total_Population_Mature_Consumers(double * Y, Parameter_Table * Table)
{
  double x;
  int i,j;

  /* Consumer subpopulation contributing to reproduction 
  */
  x = 0.0;

  assert(Table->TYPE_of_MODEL == 3); 
  
  if (Table->TYPE_of_MODEL == 3) {
    x = 0.0;
    for(j=Table->i_0; j<Table->N_E; j++){
      for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
	        if(i%Table->LOCAL_STATE_VARIABLES == (j+1)) x += Y[i]; 
    }
  }

  return(x); 
}

double Total_Population_Handling_Consumers(double * Y, Parameter_Table * Table)
{
  double x;
  int i,j;

  /* Compounds 
  */
  x = 0.0;

  assert(Table->TYPE_of_MODEL == 3); 
  
  if (Table->TYPE_of_MODEL == 3) {
    x = 0.0;
    for(j=0; j<Table->N_E; j++){
      for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
	      if( i%Table->LOCAL_STATE_VARIABLES == (j + 1 + Table->N_E) ) x += Y[i]; 
    }
  }

  return(x); 
}

double Total_Population_Triplet_Consumers(double * Y, Parameter_Table * Table)
{
  double x;
  int i,j;

  /* Triplets
  */
  x = 0.0;

  assert(Table->TYPE_of_MODEL == 3); 
  
  if (Table->TYPE_of_MODEL == 3) {
    x = 0.0;
    for(j=0; j<( Table->N_E*(Table->N_E +2 )); j++){
      for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
	      if( i%Table->LOCAL_STATE_VARIABLES == (j + 1 + 2*Table->N_E) ) x += Y[i]; 
    }
  }

  return(x); 
}

double Total_Population(double * Y, Parameter_Table * Table)
{
  double x;
  int i;

  /* This total population is calculated over all species in the system.
     It is rather a total community size if the number of species considered
     is larger than 1 
  */

  x = 0.0;
  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) x += Y[i]; 

  return(x); 
}

double Total_Population_Resource_Species (double * Y, Parameter_Table * Table)
{
  double x;
  int i, j;

  j = Table->Focal_Resource;  /*j-th state variable */
  x = 0.0;
  for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
    if(i%Table->LOCAL_STATE_VARIABLES == j) 
      x += Y[i]; 
  
  return(x);
}

double Total_Population_Species_0(double * Y, Parameter_Table * Table)
{
  double x;
  int i;

  assert(Table->TYPE_of_MODEL == 0); 
  
  if (Table->TYPE_of_MODEL == 0) {
    x = 0.0;
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++)
      if(i%2 == 0) x += Y[i]; 
  }

  return(x); 
}

double Total_Population_Species_1(double * Y, Parameter_Table * Table)
{
  double x;
  int i;

  assert(Table->TYPE_of_MODEL == 0); 
  
  if (Table->TYPE_of_MODEL == 0) {
    x = 0.0;
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) 
      if(i%2 == 1) x += Y[i]; 
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
     printf(" Type of Model is ill-defined (in Average_Individual_Density (...))\n");  
     printf(" Check input arguments. TYPE of MODEL is %d\n",
	    Table->TYPE_of_MODEL);
     printf(" The program will exit\nm");
     exit(0);
   }
   
  return (x);
}


double Total_Population_Consumers  ( double * y, Parameter_Table * Table )
{
  double x, x_0, x_1;
  int i;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
  #include <Model_Variables_Code.Include.c>
  
   if (Table->TYPE_of_MODEL == 2 || Table->TYPE_of_MODEL == 8) {
     x = 0.0;
     for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
       if(i%Table->LOCAL_STATE_VARIABLES == A)  x += y[i]; /* Free Consumers     */ 
       if(i%Table->LOCAL_STATE_VARIABLES == RA) x += y[i]; /* Handling Consumers */
     }
   }
   else if (Table->TYPE_of_MODEL == 4) {
     x = 0.0;
     for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
       if(i%Table->LOCAL_STATE_VARIABLES == A)  x += y[i]; /* Free Consumers     */ 
     }
   }
   else {
     printf(" Type of Model is ill-defined (in Total_Population_Consumers (...))\n");  
     printf(" Check input arguments. TYPE of MODEL is %d\n",
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

