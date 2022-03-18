typedef struct Master_Equationinfo
{
  
  int No_of_CONFIGURATIONAL_STATES;
  int n_x;
  int n_y;
  int n_z;
  int n_DIMENSION;

  double * Probability_Dtstribution;
  double * Probability_Distribution_Time_0; 
  double * P_n;
  double ** P_nm;
  double *** P_nml; 
  
}Master_Equation;

void Master_Equation_Allocation ( Master_Equation * ME,
				  int No_of_CONFIGURATIONAL_STATES,
				  int n_DIMENSION,
				  int n_x,
				  int n_y,
				  int n_z);

void Master_Equation_Free ( Master_Equation * ME );


void Master_Equation_Initialization (Master_Equation * ME,
				     int No_of_CONFIGURATIONAL_STATES,
				     int n_DIMENSION,
				     int n_x,
				     int n_y,
				     int n_z);

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME ); 
