#include <MODEL.h>

int Profile_Selfconsistency_Assert(int i_Strain_Sp, int * Profile, Parameter_Table * Table)
{
  int assert_true; 
  int assert_infectivity; 
  int assert_compatibility; 
  int S, k; 

    S = 0; 
    for(k = 0; k<Table->No_of_PLASMIDS; k++) 
      if( Profile[k] == 1) S++;

    if (S == 0) { /* Plasmid-free Profile */  
          assert_true = 1; 
          return(assert_true);
    }
    else {
          assert_infectivity = Infection_Condition_Assert(i_Strain_Sp, Profile, Table); 

          assert_compatibility = Compatibility_Profiles_Assert (Table, Profile, Profile, 
                                                                Table->No_of_PLASMIDS); 

          if( assert_infectivity == 1 && assert_compatibility == 1) assert_true = 1;
          else                                                      assert_true = 0; 

          return(assert_true);
    }
}
