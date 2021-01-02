#include <MODEL.h>

extern gsl_rng * r;

void Initial_Condition_from_Parameter_Table(Parameter_Table * Table, double *Y)
{
  int J;

  /* Definition of the state vector numerical order, from 0 to K, of model variables */
#include <Model_Variables_Code.Include.c>

  assert(K == Table->No_of_CELLS*Table->No_of_RESOURCES-1);

  for (J=0; J<Table->MODEL_STATE_VARIABLES; J++)
    Y[J] = Table->Vector_Model_Variables_Time_0[J];

  if(Table->RESCALING_INITIAL_TOTAL_POPULATION == 1)
    R_E_S_C_A_L_I_N_G___I_N_I_T_I_A_L___C_O_N_D_I_T_I_O_N_S ( Table, Y );

}

void R_E_S_C_A_L_I_N_G___I_N_I_T_I_A_L___C_O_N_D_I_T_I_O_N_S ( Parameter_Table * Table, double * Y)
{
    int i;
    double g_H, E_x, x, y;
    int J;

    /* Definition of the state vector numerical order, from 0 to K, of model variables */
#include <Model_Variables_Code.Include.c>

    y = 0.0;
    for( i = 0; i <= K; i++ ) y += Y[i];
    for( i = 0; i <= K; i++ ) Y[i] = Y[i]/y * Table->INITIAL_TOTAL_POPULATION;

    /* Rescaled (Corrected) Initial Conditions */
    for (J=0; J <= K; J++)
      Y[J] = Table->Vector_Model_Variables_Time_0[J] = Y[J];
}

void Initial_Condition_Centered_into_Parameter_Table (Parameter_Table * Table, double Value)
{
  /* Initial conditions from empirical data at the initial time ( -xn 0 ) */

  /* Value should represent the inital value of exposed individuals in the first age class */

  int J,n,m;
  int J_X, J_Y;
  int N_X, N_Y;

  N_X = Table->No_of_CELLS_X;
  N_Y = Table->No_of_CELLS_Y;

  for (J=0; J<Table->No_of_CELLS; J++) {
     J_X = J/N_X;
     J_Y = J%N_X;

     /* Species 0 at the center */
     if (N_X%2 == 0 && N_Y%2 == 0) {
       if ( J_X == N_X/2 && J_Y == N_Y/2 ) {
	 n=0;
	 m = Table->No_of_RESOURCES*J + n;
	 Table->Vector_Model_Variables_Time_0[m] = Value;
       }
       else {
	 n=0; 
	 m = Table->No_of_RESOURCES*J + n;
	 Table->Vector_Model_Variables_Time_0[m] = 0.0;
       }
     }
     else {

       printf(" N_X = %d\tN_Y = %d, but both N_X and N_Y must be even!!!\n",
	      N_X, N_Y);
       printf(" The program will exit\n");
       exit(0);
     }
  }

  for(n=1; n<Table->No_of_RESOURCES; n++) {
    int J_X_R = (int)(Table->No_of_CELLS_X * gsl_rng_uniform(r));
    int J_Y_R = (int)(Table->No_of_CELLS_Y * gsl_rng_uniform(r));
    
    for (J=0; J<Table->No_of_CELLS; J++) {
      J_X = J/N_X;
      J_Y = J%N_X;
      
      if ( J_X == J_X_R && J_Y == J_Y_R ) {
	
	m = Table->No_of_RESOURCES*J + n;
	Table->Vector_Model_Variables_Time_0[m] = Value;
      
      }
      else {
  
	m = Table->No_of_RESOURCES*J + n;
	Table->Vector_Model_Variables_Time_0[m] = 0.0;
      
      }
    }
  }
}

/* void Initial_Condition_Centered_into_Parameter_Table (Parameter_Table * Table, double Value) */
/* { */
/*   /\* Initial conditions from empirical data at the initial time ( -xn 0 ) *\/ */

/*   /\* Value should represent the inital value of exposed individuals in the first age class *\/ */

/*   int J,n,m; */
/*   int J_X, J_Y; */
/*   int N_X, N_Y; */

/*   N_X = Table->No_of_CELLS_X; */
/*   N_Y = Table->No_of_CELLS_Y; */

/*   for(Sp = 0; Sp < Table->No_of_RESOURCES; Sp++) { */
  
/*     for (J=0; J<Table->No_of_CELLS; J++) { */
/*       J_X = J/N_X; */
/*       J_Y = J%N_X; */
      
/*       if (N_X%2 == 0 && N_Y%2 == 0) { */
	
/* 	if ( J_X == N_X/2 && J_Y == N_Y/2 ) { */
	  
/* 	    m = Sp*J + Sp; */
/* 	    Table->Vector_Model_Variables_Time_0[m] = Value; */
	 
/* 	} */
/* 	else { */
	 
/* 	    m = Table->No_of_RESOURCES*J+n; */
/* 	    Table->Vector_Model_Variables_Time_0[m] = 0.0; */
	 
/* 	} */
/*      } */
/*       else { */
	
/* 	printf(" N_X = %d\tN_Y = %d, but both N_X and N_Y must be even!!!\n", */
/* 	       N_X, N_Y); */
/* 	printf(" The program will exit\n"); */
/*        exit(0); */
/*       } */
/*     } */
/*   } */
/* } */
