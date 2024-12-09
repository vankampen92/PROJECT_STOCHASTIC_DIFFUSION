#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void Setting_Interaction_Matrices (Parameter_Table * Table) 
{
    int i, j;

    // assert(Table->TYPE_of_MODEL == 21);  /* 21: MODEL=DIFFUSION_ECO_PLASMIDS */

    for (i=0; i<Table->No_of_STRAINS; i++) {
      
      for (j=i; j<Table->No_of_STRAINS; j++) {
        // ABB, competition matrix: symmetric competition

        if(j == i) 
          Table->ABB[i][j] = 0.0;

        else if(gsl_rng_uniform_pos(r) < Table->p_2){
          
          Table->ABB[i][j] = Table->Delta_C_0;
          Table->ABB[j][i] = Table->Delta_C_0;
        }
        
        else {
          Table->ABB[i][j] = 0.0;
          Table->ABB[j][i] = 0.0;
        } 
      }

      for (j=i; j<Table->No_of_STRAINS; j++) {
        // HBB Conjugation Matrix  
        
        if(j == i) 
          Table->HBB[i][j] = 1.0;      /* Self-conjugation between two cells of the 
                                          same bacterial type is always possible.
                                          Efective conjungation will depend on their
                                          respective plasmid profiles.  
                                          However, conjugation between two individual 
                                          cells with the same plasmid profile never 
                                          yields a transconjungant different from the 
                                          recipient. It is an event that will not 
                                          change the configuration of the system. 
                                          It will not be considered. 
                                       */
        else if (gsl_rng_uniform_pos(r) < Table->p_2){
           Table->HBB[i][j] = 1.0;
           Table->HBB[j][i] = 1.0;
        }

        else {
          Table->HBB[i][j] = 0.0;
          Table->HBB[j][i] = 0.0; 
        }
      }
         
      for (j=0; j<Table->No_of_PLASMIDS; j++) {
        // IBP Bacteria-Plasmid Interaction Matrix  (Bipartite Network)
        if(gsl_rng_uniform_pos(r) < Table->p_2)
          Table->IBP[i][j] = 1.0;
        else
          Table->IBP[i][j] = 0.0;
      }
    }

    // CPP, Plasmid-Plasmid Compatibility Matrix
    for (i=0; i<Table->No_of_PLASMIDS; i++) {

      for (j=i+1; j<Table->No_of_PLASMIDS; j++) {
        if(gsl_rng_uniform_pos(r) < Table->p_2) {
          Table->CPP[i][j] = 1.0; 
          Table->CPP[j][i] = 1.0;
        }
        else {
          Table->CPP[i][j] = 0.0;
          Table->CPP[j][i] = 0.0;
        }
      }

      Table->CPP[i][i] = 1.0; /* Reinforicing plasmid self-compatibility */
    } 
}
