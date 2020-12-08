#include <MODEL.h>

void assert_right_model_definition( Parameter_Table * P )
{

#if defined DIFFUSION
  
    assert ( P->TYPE_of_MODEL == 0 );
  
#endif

}
