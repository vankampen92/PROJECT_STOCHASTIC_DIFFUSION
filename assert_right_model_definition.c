#include <MODEL.h>

void assert_right_model_definition( Parameter_Table * P )
{

#if defined DIFFUSION
  
    assert ( P->TYPE_of_MODEL == 0 );
  
#endif
#if defined DIFFUSION_S_RESOURCES
  
    assert ( P->TYPE_of_MODEL == 1 );
  
#endif
#if defined DIFFUSION_1R1C
  
    assert ( P->TYPE_of_MODEL == 2 );
  
#endif
}
