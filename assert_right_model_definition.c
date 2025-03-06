#include <MODEL.h>

/* This is only a simple function to assert that 
 * the environmetal variable defining the model  
 * at pre-compilation time will match the number tag 
 * defined through the input argument list for every 
 * model (-y0 [MODEL NUMBER TAG])
 *
 * It also allows to give you an idea of the different 
 * models coded so far in a quick glimpse.  
 */

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
#if defined DIFFUSION_HII_2D
  
    assert ( P->TYPE_of_MODEL == 9 );
  
#endif
#if defined DIFFUSION_STOLLENBERG_3D 
  
    assert ( P->TYPE_of_MODEL == 10 );
  
#endif
#if defined DIFFUSION_HII_AC_2D
  
    assert ( P->TYPE_of_MODEL == 11 );
  
#endif
#if defined  DIFFUSION_HII_1D
  
    assert ( P->TYPE_of_MODEL == 12 );
  
#endif
#if defined DIFFUSION_BD_2D
  
    assert ( P->TYPE_of_MODEL == 13 );
  
#endif
#if defined DIFFUSION_BD_3D
  
    assert ( P->TYPE_of_MODEL == 14 );
  
#endif
#if defined DIFFUSION_STOLLENBERG_4D
  
    assert ( P->TYPE_of_MODEL == 15 );
  
#endif
#if defined  DIFFUSION_HII_nD
  
    assert ( P->TYPE_of_MODEL == 16 );
  
#endif
#if defined  DIFFUSION_AZTECA_4D
  
    assert ( P->TYPE_of_MODEL == 17 );
  
#endif
#if defined  DIFFUSION_AZTECA_4D_0
  
    assert ( P->TYPE_of_MODEL == 18 );
  
#endif
#if defined  DIFFUSION_AZTECA_4D_1
  
    assert ( P->TYPE_of_MODEL == 19 );
  
#endif
#if defined  DIFFUSION_ECOEVO_PLANTS
  
    assert ( P->TYPE_of_MODEL == 20 );
  
#endif
#if defined  DIFFUSION_ECO_PLASMIDS
  
    assert ( P->TYPE_of_MODEL == 21 );
  
#endif
#if defined  DIFFUSION_ECO_1B1P
  
    assert ( P->TYPE_of_MODEL == 22 );
  
#endif
}
