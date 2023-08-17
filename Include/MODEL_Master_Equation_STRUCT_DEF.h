typedef struct configurationinfo
{
  int * n;    /* Pointer to a configuration of the type (n_1, ..., n_D) of dimension D */
  int * anUp; /* Vector of indeces to neighgoring configurations of dimension D + 1    */
  int * anDw; 
              /* They are D+1 long to be able to store the actual relevant indeces at the
                 last position of the vector 
              */
              /* anUp[D] = a, where a is the no of UP neighboring configurations, a <= D   */
              /* anUp[D] = a, where a is the no of DOWN neighboring configurations, a <= D */
}configuration;

typedef struct Master_Equationinfo
{
  /* Up to dimension ME_n_DIMENSIONM_MAXIMUM as defined in MODEL.h */
  
  int No_of_CONFIGURATIONAL_STATES;
  int n_x;
  int n_y;
  int n_z;
  int * n_D; 
  int n_DIMENSION;

  double ** Evolving_Distribution; 
  double * Probability_Distribution;
  double * Probability_Distribution_Time_0; 
  double * P_n;
  double * PS_n;        /* Stationary Distribution */
  double ** P_nm; 
  double ** PS_nm;      /* Stationary Distribution */
  double *** P_nml;
  double *** PS_nml;    /* Stationary Distribution */
  double **** P_nmlk;
  double **** PS_nmlk;    /* Stationary Distribution */
  double ***** P_nmlkj;
  double ***** PS_nmlkj;    /* Stationary Distribution */
  double ****** P_nmlkji;
  double ****** PS_nmlkji;    /* Stationary Distribution */
  double ******* P_nmlkjih;
  double ******* PS_nmlkjih;    /* Stationary Distribution */
  double ******** P_nmlkjihg;
  double ******** PS_nmlkjihg;    /* Stationary Distribution */
  double ********* P_nmlkjihgf;
  double ********* PS_nmlkjihgf;    /* Stationary Distribution */
  double ********** P_nmlkjihgfe;
  double ********** PS_nmlkjihgfe;    /* Stationary Distribution */

  double * P_n_Marginal;
  double * P_m_Marginal;
  double * P_l_Marginal; 

  double * PS_n_Marginal;
  double * PS_m_Marginal;
  double * PS_l_Marginal; 

  /* Theoretical Marginal Probability Distributions */
  double *** MPD_T ;   /* MPD: First index represents the distribution for each dimension */
  
  /* Marginal Probability Distributions */
  double ** MPD ;   /* MPD: First index represents the distribution for each dimension */
  /* Marginal Probability Distributions at stationarity */
  double ** MPD_S;  /* MPD_S */ 
  char   ** Marginal_Probability_Label;
  
  void * Table;

  double * Vector_Model_Variables; 

  int ** TabConfi;       /* Configuration Table */
  configuration ** Co;   /* Associated Configuration Table Network Structure */

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

void Master_Equation_Configurational_State_Setup (Master_Equation * ME );

void Probability_Distribution_Vector_into_Matrix_Form( Master_Equation * ME ); 

void Master_Equation_Stationary_Distribution_Allocation ( Master_Equation * ME,
							 int No_of_CONFIGURATIONAL_STATES,
							 int n_DIMENSION,
							 int n_x,
							 int n_y,
							 int n_z );

void Master_Equation_Stationary_Distribution_Free ( Master_Equation * ME ); 
