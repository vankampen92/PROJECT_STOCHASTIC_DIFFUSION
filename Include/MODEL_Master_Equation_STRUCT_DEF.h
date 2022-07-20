typedef struct Master_Equationinfo
{
  /* Up to dimension 3 */
  
  int No_of_CONFIGURATIONAL_STATES;
  int n_x;
  int n_y;
  int n_z;
  int n_DIMENSION;

  double * Probability_Distribution;
  double * Probability_Distribution_Time_0; 
  double * P_n;
  double * PS_n;        /* Stationary Distribution */
  double ** P_nm; 
  double ** PS_nm;      /* Stationary Distribution */
  double *** P_nml;
  double *** PS_nml;    /* Stationary Distribution */

  double * P_n_Marginal;
  double * P_m_Marginal;
  double * P_l_Marginal; 

  double * PS_n_Marginal;
  double * PS_m_Marginal;
  double * PS_l_Marginal; 

  char ** Marginal_Probability_Label;
  
  void * Table;

  double * Vector_Model_Variables; 
  
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

void Master_Equation_Stationary_Distribution_Allocation ( Master_Equation * ME,
							 int No_of_CONFIGURATIONAL_STATES,
							 int n_DIMENSION,
							 int n_x,
							 int n_y,
							 int n_z );

void Master_Equation_Stationary_Distribution_Free ( Master_Equation * ME ); 
