/* This code produces 
   
   A. One single stochastic replicate of the consumer-resource model 
   in Stollenberg etat (2017)

   B. The numerical integration of the corresponding ODE system 
   
   Files are saved (nom).
   
   if CPGPLOT_REPRESENTATION is defined (in makefile or at invoking make) 
   then graphical output compare the stochastic and the deterministic 
   temporal evolution

   Compilation:

   . ~$ make

   or

   . ~$ make CPG=CPGPLOT_REPRESENTATION

   Execution (without command line arguments)
   
   . ~$ ./cotton 

   Or execution (with command line arguments)
   
   . ~$ ./cotton -pK [VALUE] -pBC [VALUE] -pBR [VALUE] -pU [VALUE] -pA [VALUE] -pDC [VALUE] -pDR [VALUE]
*/
#include <MODEL.h>

#define CURRENT_TIME_RDN_SEED
#define VERBOSE
#define No_of_EVENTS 100000

#include <global.h>

void Update_Total_Rate(int , double , double , double , double ,
		       double , double , double * , double * ,
		       int , int , int );
gsl_rng * r;

int main(int argc, char **argv)
{
  int i, i0;
  int s;

#include <default.c>

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc, argv);
  
  // #include <gsl_random_number_Setup.c>
  
  /* B E G I N : Random init */
  const gsl_rng_type * T;
  const gsl_rng_type **t, **t0;

  gsl_rng_env_setup();
  
  T = gsl_rng_taus2; //T = gsl_rng_default; T = gsl_rng_ranlux389;

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

  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); 
  /*   END: Checking Random Number Generator Setup */
  /* E N D : Init Random */
  ///////////////////////////////////////////////////////////////////

#if defined CPGPLOT_REPRESENTATION
  Parameter_CPGPLOT * CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( 1,
							      No_of_EVENTS+1,
							      0, CPG_DRIVER_NAME );  
#endif

  ///////////////////////////////////////////////////////////////////
  double time, time_0;
  double b, dx, a, v, by, dy;
  int K;

  /* This default values are not necessary any more. 
     Look at the default.c to check that default values 
     are correctly assigned to the 7 global model parameters:
     K_0, Beta_R, Delta_R, Alpha_C, Nu_C, Beta_C, Delta_C 
  
     K  = 1000;   
     b  = 10.0;   
     dx = 2.0;    
     a  = 5.0;    
     v  = 2.5;    
     by = 5.0;    
     dy = 1.0;    
  */
  
  K  = (int)K_0;     /* Carrying Capacity */
  b  = Beta_R;  /* Per Capita Growth Rate of Resources */
  dx = Delta_R; /* Per Capita Death Rate of Resources  */
  a  = Alpha_C; /* Attach Rate */
  v  = Nu_C;    /* 1/Nu_C is the Averge Handling Time    */
  by = Beta_C;  /* Per Capita Growth Rate of Consumers */
  dy = Delta_C; /* Per Capita Death Rate of Consumers  */ 
  
  /* Allocating Memmory of Parameter Structure */
  Parameter_Model * Parameters = (Parameter_Model *)calloc(1, sizeof(Parameter_Model));
  /* Initializing Param structure with parameter values */
  Parameters->K_0     = (double)K;  /* K_0     */ /* Carrying Capacity */
  Parameters->Beta_R  = b;          /* Beta_R  */ /* Per Capita Growth Rate of Resources */
  Parameters->Delta_R = dx;         /* Delta_R */ /* Per Capita Death Rate of Resources  */
  Parameters->Alpha_C = a;          /* Alpha_C */ /* Attach Rate */
  Parameters->Nu_C    = v;          /* Nu_C    */ /* 1/Nu is the Averge Handling Time    */
  Parameters->Beta_C  = by;         /* Beta_C  */ /* Per Capita Growth Rate of Consumers */
  Parameters->Delta_C = dy;         /* Delta_C */ /* Per Capita Death Rate of Consumers  */
    
  int n, m, l, n0, m0, l0; 
  
  /* Initial Conditions */
  n0 = n = 1000;
  m0 = m = 500;
  l0 = l = 0;
  time_0 = 0.0;
 
  double p0, p1, p2, p3, p4, p5, p6;
  double R;
  double * prob = (double *)calloc(7, sizeof(double) );
  
  p0 = b*(1 - n/(double)K)*(double)n;         /* Resource Growth */
  p1 = dx*(double)n;                          /* Resource Death  */
  p2 = a*(double)n /(double)K *(double)m;     /* Consumer Attack */
  p3 = v*(double)l;                           /* Consumer Handling */
  p4 = by*(double)l;                          /* Consumer Growth */
  p5 = dy*(double)m;                          /* Free Consumer Death */
  p6 = dy*(double)l;                          /* Handling Consumer Death */

  R = p0 + p1 + p2 + p3 + p4 + p5 + p6;

  prob[0] = p0/R;                     /* Resource Growth */
  prob[1] = p1/R;                     /* Resource Death  */
  prob[2] = p2/R;                     /* Consumer Attack */
  prob[3] = p3/R;                     /* Consumer Handling */
  prob[4] = p4/R;                     /* Consumer Growth */
  prob[5] = p5/R;                     /* Free Consumer Death */
  prob[6] = p6/R;                     /* Handling Consumer Death */

  /* BEGIN: Checking Discrete Sampling */
  for(i=0; i<20; i++) {
    
    s = Discrete_Sampling(prob, 7);
    
    printf("Event(%d): %d\n", i, s);
  }
  /* E N D : */
  printf("Initiating Stochastic Temporal Evolution...\n"); 
  getchar(); /* Press any key to continue... */


  /* B E G I N : Generation Stochastic Realization */ 
  printf("Initial configuration at initial time:\n");

  n = n0;
  m = m0;
  l = l0;
  time = time_0;

  printf("Initial Condition:\n");
  printf("Time = %g\t", time); 
  printf("Resource n = %d\t", n);
  printf("Serchers m = %d\t", m);
  printf("Handlers l = %d\n", l);
  printf("\n"); 

  /* Temporal_Evolution_Stochastic.dat:
     A file to save the temporal evolution of system configuration:
                     
                            (time, n, m, l)			
     by columns 
  */
  FILE * fp = fopen("Temporal_Evolution_Stochastic.dat", "w");

  double * Time_Vector = (double *)calloc(No_of_EVENTS+1, sizeof(double)); 

  fprintf(fp, "%g\t%d\t%d\t%d\n", time_0, n0, m0, l0);
  
  Time_Vector[0] = time_0;

#if defined CPGPLOT_REPRESENTATION
  double ** y_Data = (double **)calloc( 3, sizeof(double *) );
  y_Data[0] = (double *)calloc( No_of_EVENTS+1, sizeof(double) );
  y_Data[1] = (double *)calloc( No_of_EVENTS+1, sizeof(double) );
  y_Data[2] = (double *)calloc( No_of_EVENTS+1, sizeof(double) );
  
  y_Data[0][0] = (double)n0;
  y_Data[1][0] = (double)m0;
  y_Data[2][0] = (double)l0;
#endif
  
  for(i = 0; i < No_of_EVENTS; i++){
    
    time += (-1/R * log ( gsl_rng_uniform_pos(r) )); 
      
    s = Discrete_Sampling(prob, 7);

    switch(s){
    case 1 : /* Resource Growth */
      n = n + 1 ;             
      break;
    case 2 : /* Resource Death  */
      n = n - 1 ;             
      break;
    case 3 : /* Consumer Attack */
      n = n - 1;  m = m - 1;  l = l + 1;  
      break;
    case 4 : /* Consumer Handling */
      m = m + 1;  l = l - 1;
      break;
    case 5:  /* Consumer Growth */
      m = m + 1;
      break;
    case 6:  /* Free Consumer Death */
      m = m - 1;
      break;
    case 7:  /* Handling Consumer Death */
      l = l - 1;
      break;
    }
	
    printf("Time = %g\t", time); 
    printf("Resource n = %d\t", n);
    printf("Serchers m = %d\t", m);
    printf("Handlers l = %d\n", l);
    printf("\n"); 
    
    Update_Total_Rate(K, b, dx, a, v, by, dy, prob, &R,
		      n, m, l);

    fprintf(fp, "%g\t%d\t%d\t%d\n", time, n, m, l);

    Time_Vector[i+1] = time;

#if defined CPGPLOT_REPRESENTATION
    y_Data[0][i+1] = (double)n;
    y_Data[1][i+1] = (double)m;
    y_Data[2][i+1] = (double)l;
#endif

  }
  fclose(fp); 
  free(prob);
  /*     E N D : Generation Stochastic Realization */ 

#if defined CPGPLOT_REPRESENTATION
  CPG->Title[0] = '\0'; // Effectively deleting previsously defined title
  char * p = strcat(CPG->Title, "Temporal Evolution");
    
  int SAME_PLOT = 0;   
  CPG->color_Index  = 2;
  CPG->type_of_Line = 1;
  CPG->type_of_Width = 8;  
  CPG->type_of_Symbol = 1;
  /* Fixing scale of vertical (Y) axes */
  /* The X axis is fixed automatically */
  CPG->CPG_SCALE_Y   = 1;
  CPG->CPG_RANGE_Y_0 = 0.0;
  CPG->CPG_RANGE_Y_1 = 2.0 * Parameters->K_0;
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME_PLOT,
						        No_of_EVENTS+1,
							Time_Vector, y_Data[0],
							"Time",
							"Population",
							CPG->Title,
							CPG->CPG_SCALE_X,
							CPG->CPG_SCALE_Y );
  getchar();
  
  SAME_PLOT = 1;   
  CPG->color_Index  = 4;
  CPG->type_of_Line = 1;
  CPG->type_of_Width = 8;  
  CPG->type_of_Symbol = 1;  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME_PLOT,
						        No_of_EVENTS+1,
							Time_Vector, y_Data[1],
							"Time",
							"Population",
							CPG->Title,
							CPG->CPG_SCALE_X,
							CPG->CPG_SCALE_Y );
  getchar();
  
  SAME_PLOT = 1;   
  CPG->color_Index  = 5;
  CPG->type_of_Line = 1;
  CPG->type_of_Width = 8;  
  CPG->type_of_Symbol = 1;  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME_PLOT,
						        No_of_EVENTS+1,
							Time_Vector, y_Data[2],
							"Time",
							"Population",
							CPG->Title,
							CPG->CPG_SCALE_X,
							CPG->CPG_SCALE_Y );
  getchar();
#endif
  
  /* B E G I N : Numerical Integration */ 
  double * Y = (double *)calloc(3, sizeof(double)); 

  /* Initial Conditions */
  Y[0] = (double)n0;
  Y[1] = (double)m0;
  Y[2] = (double)l0;

#if defined CPGPLOT_REPRESENTATION
  y_Data[0][0] = (double)n0;
  y_Data[1][0] = (double)m0;
  y_Data[2][0] = (double)l0;
#endif
  
  fp = fopen("Temporal_Evolution_Deterministic.dat", "w"); 
  fprintf(fp, "%g\t%g\t%g\t%g\n", Time_Vector[0], Y[0], Y[1], Y[2]);

 
  double Time_Current;
  Time_Current = Time_Vector[0];
  
  for(i = 0; i < No_of_EVENTS; i++){
  
    numerical_Integration_Driver( Parameters,
				  3, Y, 
				  Time_Vector[i], Time_Vector[i+1],
				  &Time_Current );

    printf("Time: %g\t (%g, %g, %g)\n", Time_Vector[i+1], Y[0], Y[1], Y[2]);

#if defined CPGPLOT_REPRESENTATION
    y_Data[0][i+1] = Y[0];
    y_Data[1][i+1] = Y[1];
    y_Data[2][i+1] = Y[2];
#endif
    /* Comparing Time Current and Final Time for each integration interval */
    printf("Time Current (%g) vs Time_Vector[i+1] (%g)\n",
	   Time_Current, Time_Vector[i+1]);

    fprintf(fp, "%g\t%g\t%g\t%g\n", Time_Vector[i+1], Y[0], Y[1], Y[2]);
  }
  fclose(fp); 
  /*     E N D : Numerical Integration */ 

#if defined CPGPLOT_REPRESENTATION
  SAME_PLOT = 1;   
  CPG->color_Index  = 2;
  CPG->type_of_Line = 2;
  CPG->type_of_Width = 10;  
  CPG->type_of_Symbol = 1;  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME_PLOT,
						        No_of_EVENTS+1,
							Time_Vector, y_Data[0],
							"Time",
							"Population",
							CPG->Title,
							CPG->CPG_SCALE_X,
							CPG->CPG_SCALE_Y );
  getchar();
  
  SAME_PLOT = 1;   
  CPG->color_Index  = 4;
  CPG->type_of_Line = 2;
  CPG->type_of_Width = 10;  
  CPG->type_of_Symbol = 1;  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME_PLOT,
						        No_of_EVENTS+1,
							Time_Vector, y_Data[1],
							"Time",
							"Population",
							CPG->Title,
							CPG->CPG_SCALE_X,
							CPG->CPG_SCALE_Y );
  getchar();
  
  SAME_PLOT = 1;   
  CPG->color_Index  = 5;
  CPG->type_of_Line = 1;
  CPG->type_of_Width = 10;  
  CPG->type_of_Symbol = 1;  
  CPGPLOT___X_Y___P_L_O_T_T_I_N_G___S_A_M_E___P_L_O_T ( CPG, SAME_PLOT,
						        No_of_EVENTS+1,
							Time_Vector, y_Data[2],
							"Time",
							"Population",
							CPG->Title,
							CPG->CPG_SCALE_X,
							CPG->CPG_SCALE_Y );
  getchar();
#endif

#if defined CPGPLOT_REPRESENTATION
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E(CPG, 1 );
  cpgclos();
  #include <include.CPG.default.free.c>
  #include <include.FILES_to_READ.default.free.c>
  free(y_Data[0]);
  free(y_Data[1]);
  free(y_Data[2]);
  free(y_Data);
#endif
  
  free(Y);
  free(Time_Vector);
  free(Parameters);

  printf("End of program\n");
  return 0;
}

void Update_Total_Rate(int K, double b, double dx, double a, double v,
		       double by, double dy, double * prob, double * R,
		       int n, int m, int l)
{		       
  double N, M, L;

  N = (double)n;
  M = (double)m;
  L = (double)l;
  
  double p0, p1, p2, p3, p4, p5, p6;
  
  p0 = b*(1-N/(double)K)*N;    /* Resource Growth */
  p1 = dx*N;                   /* Resource Death  */
  p2 = a*(N/(double)K)*M;      /* Consumer Attack */
  p3 = v*L;                    /* Consumer Handling */
  p4 = by*L;                   /* Consumer Growth */
  p5 = dy*M;                   /* Free Consumer Death */
  p6 = dy*L;                   /* Handling Consumer Death */

  * R = p0 + p1 + p2 + p3 + p4 + p5 + p6;

  prob[0] = p0/(*R);           /* Resource Growth */
  prob[1] = p1/(*R);           /* Resource Death  */
  prob[2] = p2/(*R);           /* Consumer Attack */
  prob[3] = p3/(*R);           /* Consumer Handling */
  prob[4] = p4/(*R);           /* Consumer Growth */
  prob[5] = p5/(*R);           /* Free Consumer Death */
  prob[6] = p6/(*R);           /* Handling Consumer Death */
}
