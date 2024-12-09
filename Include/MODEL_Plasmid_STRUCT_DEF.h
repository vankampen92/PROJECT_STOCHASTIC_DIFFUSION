typedef struct Plasmidinfo
{
  int ID;               /* Identification number                 */
  int n;                /* Number of Plasmids across strains     */
  double  Cost;         /* Reproduction Cost (between 0 and 1)   */
  double  Resistance;   /* Resistance to Stress (between 0 an 1):
                           Full Resistance (1)
                           Full Sensitivity (0)
                        */
  int * Compatibility_Profile; /* List of plasmids this focal one is
                                                 compatible with!!! 
                               */    
}Plasmid;


