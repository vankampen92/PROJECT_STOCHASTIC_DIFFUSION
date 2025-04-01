#include <MODEL.h>

extern gsl_rng * r;   /* Global generator (define at the main program level */
#define RANDOM gsl_rng_uniform_pos(r)

/* 
   This function advances the system one single time step, by changing
   system configuration by one single discrete event.  When 
   the stochastic dynamics is implemented following Gillespie Direct 
   algorithm the exact time at which this occurs is sampled from an exponential
   distribution at a rate that represents the total rate at which the
   system changes configuration at any given time (Stochastic_Rate->Total_Rate).
   When it is implemented using the next reaction method, time increases
   according to the shortest waiting time over all events. These times are 
   sampled from the specific expential distribution associated to each event.  

   This function calls two fundamental functions: Exectute_One_Step() and 
   either Temporal_Dynamics() or Temporal_Dynamics_Update(), where the  
   model-specific stochastic dynamics and rates are defined.
*/

int Advance_Current_Time( Parameter_Table * Table, 
			                    Stochastic_Rate * Rate, double * Time_Current, 
                          int * New )
{
  /* (*New) counts the number of infection events */ 
  int k;
  double inter_event_time; 
  int Event; 
  double Max_Probability     = Rate->max_Probability;     
  int no_Patch               = Table->No_of_CELLS;
  /* Patch[0] and Patch[1] recover two patch labels of either: 
     1. Patch[0] suffers from outmigration and Patch[1] received an immigrated individual 
     or 
     2. When Patch[0] is equal to Patch[1], the event occurred in this patch did not 
     involve any movement between patches
  */
  
  int * Patch                = (int *)calloc(5, sizeof(int));
  Patch[0]                   = 0; /* Single central cell system (always true when No_of_CELLS = 1 (-HM 1) */
  Patch[1]                   = 0; /* Single central cell system (always true when No_of_CELLS = 1 (-HM 1) */
  Patch[2]                   = Table->No_of_RESOURCES; /* Impossible Species ID */
  Patch[3]                   = Table->No_of_RESOURCES; /* Impossible Species ID */
  /* (Sp_ID equal to Table->No_of_RESOURCES is impossible, because 0 <= Sp_ID <= Table->No_of_RESOURCES-1)*/
  Patch[4]                   = 1;                      /* 0: No configurational change has occured */
                                                       /* 1: A configuration change has occured    */
                                                       
  /* This array stores information about the event that has occurred: 
     Patch[0] Local Patch (where the event has occurred ) Default: 0, when No_of_CELLS = 1 (-HM 1) 
     Patch[1] Patch receiving the immigrant if a movement event has occurred. Default: 0, 0, when No_of_CELLS = 1 (-HM 1) 
     Patch[2] Index of an extra species involved in the event
     Patch[3] Index of a second extra species also involved in the event
     Patch[4] Flag indicating that a configurational change has occured. 
             0: No configurational change has occured
             1: A configuration change has occured, and the update of the array of rates is required
  */
  
  Parameter_Model * P        = Table->P; 
  Community ** Village       = Table->Patch_System;   

  int * J                    = Table->Vector_Model_Int_Variables;
  double * Y                 = Table->Vector_Model_Variables;
  
  int DISCARTING_EXTINCTIONS = Table->T->DISCARTING_EXTINCTIONS;  
  double time_Factor         = Table->T->Time_Scale_Unit;

  if( (Rate->Total_Rate) > 0){

    #ifndef PRIORITY_QUEU_SUPER_OPTIMIZATION
    #ifndef REUSE_RANDOM_NUMBER
      inter_event_time = -1./(Rate->Total_Rate) * log(RANDOM);
      (*Time_Current) += inter_event_time;
    #endif
    #endif 
    /* BEGIN : Stochastic Dynamic is actually performed : 
               The Community "Village" is updated accordingly 
    */
    Execute_One_Step( Village, Table, Max_Probability, &Event, Patch );
    /*   END : Stochasctic Dynamics * * * * * */

    // if(Event == 0) (*New)++; /* Accumulating events of Type 0 between times */
    #if defined REUSE_RANDOM_NUMBER
      #if defined BINARY_TREE_SUPER_OPTIMIZATION
        inter_event_time = -1./(Rate->Total_Rate) * log(Rate->Reusable_Random_Number);
        (*Time_Current) += inter_event_time;
      #else
        printf(" The reuse of random number is only possible if the stochastic opitmization\n");
        printf(" at work is BINARY_TREE_SUPER_OPTIMIZATION, but this optimization level has\n");
        printf(" not been defined at compilation time (through the makefile)\n");
        printf(" The program will exit.");
        exit(0);
      #endif
    #endif

    #if defined PRIORITY_QUEU_SUPER_OPTIMIZATION
      /* Advance time... */
      /* Tree_Node_Index: Priority array of indexed pointers to all tree nodes */
      int Index_Node_x = Table->TOTAL_No_of_EVENTS * Patch[0] + Event;
      
      if( Table->Tree_Node_Index[Index_Node_x] != Table->Treeroot ) {
	     printf("Table->TOTAL_No_of_EVENTS = %d\n", Table->TOTAL_No_of_EVENTS); 
	     printf("Event    = %d\n", Event); 
	     printf("Patch[0] = %d\n", Patch[0]);
      }

      assert(Table->Tree_Node_Index[Index_Node_x] == Table->Treeroot);
      
      (*Time_Current)  = Table->Tree_Node_Index[Index_Node_x]->value;
      
      Rate->Stochastic_Time = (*Time_Current);   
    /* --------------- */
    #endif 
  }
  else{    

    printf(" Total System Probability Rate of Change = %g\n",  Rate->Total_Rate);
    printf(" When the total Rate is zero, the average time to the next event is infinite!!!\n");
    printf(" So, events can no longer occur. The system does not react any more.\n");
    
    return(1);  
  }
  
  /* BEGIN: Updating of the Total Rate of Change for the system       */
  if (Table->T->TYPE_of_TIME_DEPENDENCE > 0) {
    /* Update_Time_Dependence (Table->TYPE_of_TIME_DEPENDENCE, Time_Current, Table ); */
    // Trend_Time_Dependence( Table, (*Time_Current) );
  }
  
  #if defined STOCHASTIC_OPTIMIZATION 
    Temporal_Dynamics_Update( Village, Table, Rate, Event, Patch );
  #else
    Temporal_Dynamics(Village, Table, Rate);
  #endif


  #if defined VERBOSE
  if( * Time_Current < 1.0) { 
    /* Only initial times are printed out...just for a check!!! */
    printf("Time = %g\t Type of Event = %d in Patches (%d, %d)\n", 
            * Time_Current, Event, Patch[0], Patch[1]);
    if (Table->No_of_CELLS == 1) 
      Print_Meta_Community_Patch_System (Table);
      
    printf("\n");
  } 
  #endif

  #if defined VERBOSE
  if( * Time_Current >= 1.0) { 
    printf("Time = %g\t Type of Event = %d in Patches (%d, %d)\n", 
          * Time_Current, Event, Patch[0], Patch[1]);
    printf(" Current system configurtion (t = %g):\n",  (*Time_Current) );
    Print_Meta_Community_Patch_System (Table);
    printf("\n");
    getchar();
  }
  #endif

  free(Patch);
  return(0);
}
