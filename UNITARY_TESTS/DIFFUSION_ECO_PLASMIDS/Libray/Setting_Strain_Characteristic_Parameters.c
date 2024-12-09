#include <MODEL.h>

void Setting_Strain_Characteristic_Parameters (Parameter_Table * Table)
{
    int i, j, k, n; 
    double c, COST, RESISTANCE;  

    n = 0;
    for (i=0; i<Table->No_of_STRAINS; i++) {  
      for(j=0; j<Table->n[i]; j++) { 
        
        c = 0.0; COST = 1.0; RESISTANCE = 0.0;        
        for(k=0; k<Table->No_of_PLASMIDS; k++) {
          if(Table->Strain_Profiles[i][j][k] == 1) { 
            c = Table->Alpha_C[k];
            RESISTANCE = MAX(RESISTANCE, Table->Nu_C[k]);
            COST *= (1.0 - c);
          }
        }

        Table->Beta_AP[n]  = Table->Beta_R * COST;                                        /* Bacteria Cell Division Rates */
        Table->Eta_RP[n]   = Table->Lambda_R_1;                                           /* Bacteria Conjugation Rates   */
        Table->Delta_AP[n] = Table->Delta_R_0 + Table->Delta_R_1 * (1.0 - RESISTANCE);    /* Bacteria Death Rates         */
        Table->Mu_RP[n]    = Table->Mu;                                                   /* Bacteria Diffusion Rates     */
        Table->Segregation_Error[n]= Table->p_1;                                          /* Bacterial Segregation Error  */

        n++; 
      }
    }

    assert(n == Table->No_of_RESOURCES);
}  
