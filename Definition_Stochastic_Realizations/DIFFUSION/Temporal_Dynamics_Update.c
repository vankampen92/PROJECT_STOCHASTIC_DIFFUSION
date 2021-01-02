/* Important: Under Construction: Check updating, please!!!!!!!!!!!!!!!!!!!  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2010 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

void Temporal_Dynamics_Update( Community ** My_Community,
			       Parameter_Table * Table, Stochastic_Rate * Rate,
			       int Type_of_Event, int * Patch)
{
  /* This function calculates the stochastic rates after the execution of a stochastic event
     in terms of the old ones, with no recalculation. This is a way to optimize the algorithm.
     It is always worth trying this optimization in very sparsely coupled systems. 
  */
  /* Input arguments:
     
     . My_Community is a pointer to the whole patch system
     . Table        is a pointer to a Parameter_Table type of structure from which other 
                    structures hang. 
     . Rate         is a pointer to the structure Stocastic Rate 
     . Type_of_Event is a label to the event occurring in a Patch 
     . Patch        is an array containing the two patches involved in a movement event. 
                    Patch[0] is the patch sending the individual
		    Patch[1] is the patch receiving the individual. 

     Output arguments: 
     . Rate         Stochastic Rate is updated from previous value (without recalculating)
  */
  
  int x, y; 
  int Sp; 
  Parameter_Model * P = Table->P; 
  Community * Pa;
  double OutMigration; 
  
  x = Patch[0]; y = Patch[1];
  Sp = Type_of_Event%Table->No_of_RESOURCES;
  
  /* Patch sending one individual out (x ---> y) */
  Pa    = My_Community[x];

  OutMigration = Pa->Total_Per_Capita_Out_Migration_Rate[Sp];

  assert( 4 * Table->Mu ==  OutMigration ); 
        
  Pa->rate[Sp] = OutMigration;   Pa->rToI[Sp]  -= OutMigration;
    
                                 Pa->ratePatch -= OutMigration;

                              Rate->Total_Rate -= OutMigration;
			      
  Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch );
  
  /* Patch receiving one individual in (x ---> y) */
  Pa    = My_Community[y];

  OutMigration = Pa->Total_Per_Capita_Out_Migration_Rate[Sp];

  assert( 4 * Table->Mu ==  OutMigration ); 
        
  Pa->rate[Sp] = OutMigration;   Pa->rToI[Sp]  += OutMigration;
    
                                 Pa->ratePatch += OutMigration; 

                              Rate->Total_Rate += OutMigration;
			      
  Rate->max_Probability = MAX( Rate->max_Probability, Pa->ratePatch );
 
  if(Rate->Total_Rate <= 0.0){
      printf("\n");
      printf(" R is the total temporal rate of system configuration change\n");
      printf(" R = %g\n", Rate->Total_Rate );
      printf(" As R is zero or negative, no change is possible\n");
      printf(" R shouldn't be negative. If it is, there are definitely some errors in the code\n");
      printf("\n");
      if( Rate->Total_Rate < 0.0) exit(0);
  }
}
