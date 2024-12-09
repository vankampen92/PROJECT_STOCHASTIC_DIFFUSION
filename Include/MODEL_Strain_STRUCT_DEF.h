typedef struct Straininfo
{
  /* A strain is a particular subpopulation with a given plasmid profile 
  */

  int ID_0;           /* Potential ID according to Profile, from 1 to (2^n -1)   */
  int ID;             /* True ID Realized and Reorded                            */
  int * Profile;
  int n;              /* Number of Individuals in this Strain */
  double  Death_0;    /* Per Capita (density-independent) Death Rate */
  double  Death_1;    /* Per Capita (density-independent) stress induced Death Rate */
  double  Beta;       /* Per Capita Growth Rate                      */
  double  Segregation_Error; /* Segregation Error                    */
  double  Gamma;      /* Conjugation/Encounter Rate                  */
  double  Mu_0;       /* Strain Diffusion Rate                       */

  int * Competition_List;  /* List of strains (by True ID) that cause competition-induced extra mortality 
                              on the given focal one */

  int * Conjugation_List;  /* List of strains (by True ID) a given focal strain can conjugate with */

  int * Recipient_List;    /* List of train and profiles (supbpopualtions, by true ID) that potentially 
                              may act as recipients in a conjugation event 
                           */
  int * Donor_List;         /* List of donors (by true ID) a given Strain can potentially receive plasmids 
                              from 
                           */                                      
}Strain;


