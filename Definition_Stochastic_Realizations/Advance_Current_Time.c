#include <MODEL.h>

extern gsl_rng * r;   /* Global generator (define at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

/* This function advances the system one single time step, by changing
   system configuration by one single discrete event.  Since
   the stochastic dynamics is implemented following Gillespie algorithm
   the exact time at which this occurs is sampled from an exponential
   distribution at rate the represents the total rate at which the
   system changes configuration at any given time (Stochastic_Rate->Total_Rate).

   This function calls two fundamental functions, step() and Temporal_Dynamics(), 
   where the actual model specific stochastic dynamics is defined.
*/

int Advance_Current_Time( Parameter_Table * Table, 
			  Stochastic_Rate * Rate, double * Time_Current, int * New )
{
  /* (*New) counts the number of infection events */ 
  int k;
  double inter_event_time; 
  int Event; 
  double Max_Probability     = Rate->max_Probability;     
  int no_Patch               = Table->No_of_CELLS;
  int * Patch                = (int *)calloc(2, sizeof(int)); 
  
  Parameter_Model * P        = Table->P; 
  Community ** Village       = Table->Patch_System;   

  int * J                    = Table->Vector_Model_Int_Variables;
  double * Y                 = Table->Vector_Model_Variables;
  
  int DISCARTING_EXTINCTIONS = Table->T->DISCARTING_EXTINCTIONS;  
  double time_Factor         = Table->T->Time_Scale_Unit;

  if( (Rate->Total_Rate) > 0){

    inter_event_time = -1./(Rate->Total_Rate) * log(RANDOM);
    (*Time_Current) += inter_event_time;
    
    /* BEGIN : Stochastic Dynamic is actually performed : Village is updated accoundingly */
    Execute_One_Step( Village, Table, Max_Probability, &Event, Patch );
    /*   END : Stochasctic Dynamics * * * * * */

    if(Event == 0) (*New)++; /* Accumulating events of Type 0 between times */

  }
  else{    

    printf(" Total System Probability Rate of Change = %g\n",  Rate->Total_Rate);
    printf(" When the total Rate is zero, the average time to the next event is infinite!!!\n");
    printf(" So, events can no longer occur. The system does not react any more.\n");
    
    return(1);  
  }
  
  /* BEGIN: Calculation of the Total Rate of Change for the system       */
  if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) {
    /* Update_Time_Dependence (Table->TYPE_of_TIME_DEPENDENCE, Time_Current, Table ); */
    // Trend_Time_Dependence( Table, (*Time_Current) );
  }

  // Temporal_Dynamics(Village, Table, Rate); 
  Temporal_Dynamics_Update( Village, Table, Rate, Event, Patch );

  free(Patch); 
  /*   END: Calculation of Total Rate of Change */

#if defined VERBOSE
  printf(" Total population across the system at current time (t = %g)\n", 
	 (*Time_Current) );
  
  Print_Meta_Community_Patch_System (Table);
#endif

  
  return(0);
}
  
