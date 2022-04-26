#include <MODEL.h>

#define FEEDING_EVENT_FACTOR 100

void Model_Parameters_Master_Equation(Parameter_Table * Table,
				      int * No_of_CONFIGURATIONAL_STATES,
				      int * n_DIMENSION,
				      int * n_x, int * n_y, int * n_z)
{
  double Time     = Table->T->Time_1 - Table->T->Time_0; 
  double N        = Average_Number_of_Feeding_Events( Table, Time ); 
  
  * n_DIMENSION = 2;
  * n_z         = 1;
  * n_y         = Table->TOTAL_No_of_CONSUMERS;
  * n_x         = FEEDING_EVENT_FACTOR * (int)N; 
  
  * No_of_CONFIGURATIONAL_STATES = Table->TOTAL_No_of_CONSUMERS * ( * n_x ) ;
}

double Average_Number_of_Feeding_Events( Parameter_Table * Table, double T )
{
  double n_T;
  double n_R, m_0, K_R, Theta, Nu;

  K_R   = (double)Table->K_R;
  n_R   = (double)Table->TOTAL_No_of_RESOURCES; /* which determines resource density */
  m_0   = (double)Table->TOTAL_No_of_CONSUMERS;

  Theta = Table->Alpha_C_0 * n_R/K_R;  /* -H9  [Alpha_C_0] */
  Nu    = Table->Nu_C_0;               /* -H10 [Nu_C_0]    */
  
  n_T = Theta * Nu/(Nu+Theta) * ( T + Theta/Nu/( Nu + Theta) * (1.0-exp(-(Theta+Nu)*T)) ) * m_0;
  
  return (n_T);
}

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME )
{
  int i;
  int n,m,l;
  int s; 
  /* Notice that this only works for probability distributions defined in maximum 
     three dimensions where each dimension goes from 0 to a maximum constant value: 

     n = 0, ..., ME->n_x
     m = 0, ..., ME->n_y
     l = 0, ..., ME->n_z
  */
  
  double * y = ME->Probability_Distribution; 

  if (ME->n_DIMENSION == 1) {
    for(i=0; i<ME->n_x; i++) ME->P_n[i] = y[i];
  }  
  else if (ME->n_DIMENSION == 2) {
    
    for(i=0; i<ME->No_of_CONFIGURATION_STATES; i++) {
      n = i/ME->n_y;
      m = i%ME->n_y;
      ME->P_nm[n][m] = y[i];
    }
    
  }
  else if (ME->n_DIMENSION == 3) {

    for(i=0; i<ME->No_of_CONFIGURATION_STATES; i++) {
      n = i/(ME->n_y*ME->n_z);
      s = i%(ME->n_y*ME->n_z);
      m = s/ME->n_z;
      l = s%ME->n_z
	
      ME->P_nml[n][m][l] = y[i];
    }
    
  }
  else {
    printf(" This structure is only prepared to accept maximum three dimensions, but\n");
    printf(" you probability distribution seems to have %d dimensions!!!\n", ME->n_DIMENSION); 
    assert(ME->n_DIMENSION < 4);
  }
}








