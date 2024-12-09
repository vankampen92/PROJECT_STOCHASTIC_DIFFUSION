#include <MODEL.h>

void Setting_Plasmid_Characteristic_Parameters (Parameter_Table * Table)
{
    //Plasmid Cost and Resistance 
    int i; 

    for (i=0; i<Table->No_of_PLASMIDS; i++) {
      Table->Alpha_C[i]  = Table->Alpha_C_0;    /* Plasmid reproduction costs */
      Table->Nu_C[i]     = Table->Nu_C_0;       /* Plasmid resistance         */
    }
} 
