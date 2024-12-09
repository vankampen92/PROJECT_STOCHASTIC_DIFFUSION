#include <MODEL.h> 
#include "Parameter_Definitions.h"

#define VERBOSE
#define CURRENT_TIME_RDN_SEED

#include "global.h"
gsl_rng * r; /* Global generator defined in main.c */

void ArgumentControl(int argc, char **argv);

int main(int argc, char **argv)
{
  Parameter_CPGPLOT * CPG;
  int colors[3] = {2,12,8};
  
  #include "default.c"
 /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  
  int i, j, k, l, n, N;
   
  /* GNU Scientific Library:  Random numbers set up */
  const gsl_rng_type * T;
  const gsl_rng_type **t, **t0;

  gsl_rng_env_setup();  /* Environment variables  
                          GSL_RNG_TYPE and GSL_RNG_SEED are read  
                          to set the corresponding library variables 
                          gsl_rng_default and gsl_rng_default_seed. 
                        */
  T = gsl_rng_taus2;   //T = gsl_rng_default; T = gsl_rng_ranlux389;
  r = gsl_rng_alloc(T);
  t0 = gsl_rng_types_setup ();

#if defined VERBOSE
  printf("Random number information... Available generators:\n");
  for(t=t0; *t != 0; t++){
    printf("%s, ", (*t)->name);
  }
#endif
  printf("\n In this example, the random generator at work... is the '%s' generator\n\n",
         gsl_rng_name(r));

#if defined CURRENT_TIME_RDN_SEED
  GSL_Init_Random_Seed(r);
#else
  GSL_Init_Random_Seed_from_File(r);
#endif

/* End of Random number setting */
#if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); //getchar();
  /*   END: Checking Random Number Generator Setup */
#endif

  /* Allocating and Setting Matrices:
            1.  Competition Matrix, ABB 
            2.  Conjugation Matrix, HBB 
            3.  Plasmid-Bacteria Infection Matrix, IPB
            4.  Plasmid-Plasmid Compatibility Matrix, CPP
  */
  Parameter_Table * Table = (Parameter_Table *)calloc(1, sizeof(Parameter_Table));
    
  Preparing_Initial_System_Configuration( Table );  /* Allocating and initializing data structures */

  Time_Control * Time = (Time_Control *)calloc(1, sizeof(Time_Control));  
  Time->Time_Vector   = (double *)calloc(No_of_TIMES, sizeof(double));
  Time_Control_Upload ( Time, Table );

  int No_of_PANELS = Table->No_of_RESOURCES;
  int No_of_POINTS = No_of_TIMES;
  CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( No_of_PANELS, No_of_POINTS, 0,
                                          CPG_DRIVER_NAME );
  printf(" Parameter_CPGPLOT data structure has been allocated and initiated\n");

  Table->CPG = CPG;

  int STATE = Deterministic_Time_Dynamics( Table );

  De_Allocating_Initial_System_Configuration ( Table );

  printf("\n ------------------------------------------------------------------------------ \n");
  printf("No of RESOURCES (No of STRAIN IDs): %d\n", Table->No_of_RESOURCES);
  printf("No of STRAIN TYPES: %d\n", Table->No_of_STRAINS);
  printf("No of PLASMIDS: %d\n", Table->No_of_PLASMIDS);
  printf("No of TOTAL EVENTS: %d\n", Table->No_of_STRAINS * 6 + Table->No_of_CONJUGATION_EVENTS);
  printf("No of Conjugation Pairs (Donor, Recipient): %d\n", Table->No_of_CONJUGATION_EVENTS);

  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( CPG, No_of_PANELS );
  cpgclos();

  free(Time->Time_Vector);
  free(Time);
  free(Table);     


  printf(" ---> End of Program\n");
  return(0);
}

