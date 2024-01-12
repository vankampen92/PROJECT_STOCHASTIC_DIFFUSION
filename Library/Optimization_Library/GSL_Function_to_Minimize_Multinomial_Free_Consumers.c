#include <MODEL.h>

double GSL_Function_to_Minimize_Multinomial_Free_Consumers( const gsl_vector * x, void * Par )
{
  /* This GSL function allows the evaluation of a likelihood function to minimize or 
     evaluate. 

     This likelihood corresponds to the probability of observing a number of 
     consumers feeding on multiple resource types at different times for a feeding 
     dynamics represented by Holling Type II model. 

     This is the reason there is an assert to verify the corret model has been chosen 
     and compiled. 

     Output value:
     . Value, the total multiplicative likelihood: Probability of data given a choice 
       of parameter values: 
                            -log( P(Data | parameters) ) 
   */
  int i,j;
  int n_H; 
  double Value, Theta_Sum, p_Sum, Nu_Sum; 
  double t; 

  Parameter_Fitting * F   = (Parameter_Fitting *)Par;

  int No_of_POINTS         = F->Data->No_of_POINTS;
  int No_of_VARIABLES      = F->Data->No_of_VARIABLES;
  double ** Data           = F->Data->N;
  double * Theory;
  Parameter_Table * Table  = F->Table; 
  Parameter_Space * Space  = F->Space;
  int No_of_PARAMETERS     = F->Space->No_of_PARAMETERS;
  int x_is_BOUNDED; 
                                                                      
  assert ( Table->TYPE_of_MODEL == 16);                               /* 16: DIFFUSION_HII_nD */
  assert ( Table->SUB_OUTPUT_VARIABLES == Table->No_of_RESOURCES );   /* No of RESOUCE TYPES  */
  assert ( No_of_VARIABLES == Table->No_of_RESOURCES );                                   

  Theory = (double *)calloc( No_of_POINTS, sizeof(double) ); 
    
  x_is_BOUNDED = 1; /* By default, the neg Log Likelihood is evaluated */
  
  if(F->Minimization == 1)
    x_is_BOUNDED = Checking_for_Parameter_Boundaries( F, x );
  
  if( x_is_BOUNDED == 1 ) { 
						    
    if(No_of_PARAMETERS > 0) 
      Vector_Entries_into_Parameter_Table ( x, Table,
					    Space->Parameter_Index, No_of_PARAMETERS );
    
    if(Table->No_of_IC > 0) 
      Vector_Entries_into_Parameter_Table_Initial_Condition ( x, Table, 
							      Table->IC_Space->Parameter_Index,
							      No_of_PARAMETERS,
							      Table->No_of_IC );
    if(Table->No_of_ERROR_PARAMETERS > 0)
      Vector_Entries_into_Parameter_Table_Error_Model ( x, Table,
							Table->E_Space->Parameter_Index,
							No_of_PARAMETERS,
							Table->No_of_IC,
							Table->No_of_ERROR_PARAMETERS );

    Resetting_Alpha_Nu_Vectors_Constant (Table); /* Resseing Alpha's and Nu's        */
    Resetting_Multiresource_Levels (Table);      /* Rebuilding Theta_Consumers[] !!! */

    // Writing_Alpha_Nu_Theta_Vectors(Table);
     
    int * n = (int *)calloc(No_of_VARIABLES+1, sizeof(int));
    double * p = (double *)calloc(No_of_VARIABLES+1, sizeof(double)); 

    // Theta_Nu_Sums (Table->Theta_Consumers, Table->Nu_Consumers, 
    //                &Theta_Sum, &Nu_Sum, Table->No_of_RESOURCES); 
                   
    Theta_Sum = 0.0; Nu_Sum  = 0.0;  
    for (j=0; j<Table->No_of_RESOURCES; j++){ 
      Theta_Sum += Table->Theta_Consumers[j];
      // printf("Nu[%d] = %g [Nu_C_0 = %g]  ", j, Table->Nu_Consumers[j], Table->Nu_C_0);

      assert(Table->Nu_Consumers[j] == Table->Nu_C_0); /* Analytic Formula only valid 
                                                          under this restriction 
                                                        */
      Nu_Sum += (Table->Theta_Consumers[j]/Table->Nu_Consumers[j]);
    }   
    // printf("\n");

    for(j=0; j<No_of_POINTS; j++) { 

      t = F->Data->Time_Data_Vector[j];
      
      n_H = 0; p_Sum = 0.0; 
      for(i=0; i<Table->SUB_OUTPUT_VARIABLES; i++) {
          n[i] = (int)Data[i][j];
          n_H += n[i]; 
	
        //  p[i] = Function___p_of_t (Table->Theta_Consumers[i], 
        //                            Table->Nu_Consumers[i], 
        //                            Nu_Sum, Theta_Sum, t);

	      p[i] = (Table->Theta_Consumers[i]/Table->Nu_Consumers[i])/(1.0 + Nu_Sum) * (1.0-exp(-(Table->Nu_Consumers[i] + Theta_Sum)*t));

	      p_Sum += p[i];
      }
      
      n[No_of_VARIABLES] = Table->TOTAL_No_of_CONSUMERS - n_H; /*-HN 20 [TOTAL_No_of_CONSUMERS]*/ 
      p[No_of_VARIABLES] = 1.0 -p_Sum; 

      assert(n[No_of_VARIABLES] <= Table->TOTAL_No_of_CONSUMERS);
      assert(n[No_of_VARIABLES] >= 0);

      Theory[j] =  -gsl_ran_multinomial_lnpdf(No_of_VARIABLES+1, p, n);                
    }

    free(n); free(p);

    int Theory_is_NOT_a_NUMBER = 0;
    Theory_is_NOT_a_NUMBER = da_vector_isnan(Theory, No_of_POINTS);
    if( Theory_is_NOT_a_NUMBER == 0 ) { 
      Value = 0.0;
	    for(j=0; j<No_of_POINTS; j++) Value += Theory[j]; 
    }
    else Value = DBL_MAX; /* Likelihood takes a non-feasible value */   

  }
  else Value = DBL_MAX;   /* Parameter combination outside boundaries */

  free(Theory);
  
  return(Value);
}

void Theta_Nu_Sums (double * Theta_Vector, double * Nu_Vector, 
                    double * Theta_Sum, double * Sum_ThetaNu_Ratio, 
                    int N)
{
    int j; 

    * Theta_Sum = 0.0; * Sum_ThetaNu_Ratio = 0.0;  
    for (j=0; j<N; j++){ 
        * Theta_Sum += Theta_Vector[j];
        * Sum_ThetaNu_Ratio += (Theta_Vector[j]/Nu_Vector[j]);
    }
}

double Function___poft (double Theta, double Nu, 
                        double Sum_ThetaNu_Ratio, double Theta_Sum, 
                        double t)
{
    double p; 

    p = (Theta/Nu)/(1.0 + Sum_ThetaNu_Ratio) * (1.0-exp(-(Nu + Theta_Sum)*t));

    return p;
}
