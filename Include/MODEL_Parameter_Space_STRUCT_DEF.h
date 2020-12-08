typedef struct Parameter_Spaceinfo
{
  /* Model Parameters defining the boundary of
     the Parameter Space */
  Parameter_Model * Parameter_min;
  Parameter_Model * Parameter_MAX;

  /* Model Parameters defining parameter dependent accuracies */
  Parameter_Model * Parameter_Acc;

  /* * * Model Parameters for local community dynamics  * * */
  gsl_vector * P_MAX;         /* Max and Min for each parameter */
  gsl_vector * P_min;
  gsl_vector * Accuracy;      /* Different accuracy for each searcheable parameter */


  int * Parameter_Index;

  /* This will define GRID on which different parameter combinations will be
     evalutated */
  int * N;
  double ** N_Par_Value;

  #include <include.Parameter_Space.global.h>

#if defined PARAMETRIC_CONFIGURATIONS
  int No_Parametric_Configurations;
  double ** Parametric_Configuration;
#endif
  
}Parameter_Space;

#if defined PARAMETRIC_CONFIGURATIONS
void Parameter_Space_Alloc( Parameter_Space * , int , int , int *  );
#else
void Parameter_Space_Alloc( Parameter_Space * , int , int *  );
#endif

void Parameter_Space_Initialization(Parameter_Space * , int , double, int, 
				    int * , int * , double ** , double * );
void Parameter_Space_Free(Parameter_Space *, int );


void Parameter_Index_Checking_Ordering (int * , int * , int ,
					int * , int * , int * , int * ,
					int * , int * , int * );







