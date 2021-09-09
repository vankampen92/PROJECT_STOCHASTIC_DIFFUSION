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
#if defined DIFFUSION_1RnC_E
  
    assert ( P->TYPE_of_MODEL == 3 );
  
#endif
#if defined DIFFUSION_1R1C_2D
  
    assert ( P->TYPE_of_MODEL == 4 );
  
#endif
#if defined DIFFUSION_DRAG
  
    assert ( P->TYPE_of_MODEL == 5 );
  
#endif
#if defined DIFFUSION_VRG
  
    assert ( P->TYPE_of_MODEL == 6 );
  
#endif
#if defined DIFFUSION_MR
  
    assert ( P->TYPE_of_MODEL == 7 );
  
#endif
#if defined DIFFUSION_1R1C_2D_STO_4D
  
    assert ( P->TYPE_of_MODEL == 8 );
  
#endif
}
