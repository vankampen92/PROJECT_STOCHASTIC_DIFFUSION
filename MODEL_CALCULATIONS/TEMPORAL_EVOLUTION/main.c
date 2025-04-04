/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2018 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <Include/MODEL.h>

#include "global.h"

gsl_rng * r; /* Global generator defined in main.c */

/* This code calculates the stochastic and determinisitic temporal evolution of several resource-consumer models 

   Compilation (see makefile environment variable 'MODEL'). Some compilation commands as example:

   . ~$ make
   . ~$ make STATIONARITY=STATIONARY_POINT_REPRESENTATION MODEL=DIFFUSION_STOLLENBERG_4D
   . ~$ make STATIONARITY=NON_STATIONARY_POINT_REPRESENTATION MODEL=DIFFUSION_STOLLENBERG_4D

   Exectution:
   
   1 species examples:                                       (OUTPUT_VARIABLES_GENUINE will be 4) 
   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 1 -HM 100 -HX 10 -HY 10 -Hu 0.1 \ 
                   -n 2 -v0 0 -v1 60 -G0 1 -G1 2 \ 
                   -tn 50 -t0 0.0 -t1 50.0 -t4 0 -tR 10 -xn 0 -xN 1000 -HN 1000 \
                   -G2 1 -G3 0.0 -G4 50.0 -G5 1 -G6 0.0 -G7 6000.0 

   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 -Hu 0.1 \ 
                   -n 1 -v0 5054 -G0 1 -G1 1 \
                   -tn 100 -t0 0.0 -t1 200.0 -t4 0 -tR 10 -xn 0 -xN 1000 -HN 1000 \ 
                   -G2 1 -G3 0.0 -G4 50.0 -G5 1 -G6 0.0 -G7 1100.0

   2 species examples:                                       (OUTPUT_VARIABLES_GENUINE will be 5) 
   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 2 -HM 10000 -HX 100 -HY 100 -Hu 0.5 \
                   -n 1 -v0 10105 -G0 1 -G1 1 \
                   -tn 100 -t0 0.0 -t1 30.0 -t4 0 -tR 5 -xn 0 -xN 1000 -HN 1000 \
                   -G2 1 -G3 0.0 -G4 30.0 -G5 1 -G6 0.0 -G7 1100.0

   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 2 -HM 100 -HX 10 -HY 10 -Hu 0.5 \
                   -n 1 -v0 115 -G0 1 -G1 1 \
                   -tn 20 -t0 0.0 -t1 5.0 -t4 0 -tR 5 -xn 0 -xN 98 -HN 98 \
                   -G2 1 -G3 0.0 -G4 5.0 -G5 1 -G6 0.0 -G7 100.0

   3 species example:                                        (OUTPUT_VARIABLES_GENUINE will be 6) 
   .~$ ./DIFFUSION -y0 0 -y2 1 -HS 3 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -n 1 -v0 15156 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 30.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 30.0 -G5 1 -G6 0.0 -G7 1100.0

   1 species (with external immigration and death) example:  (OUTPUT_VARIABLES_GENUINE will be 4) 
   .~$ ../DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.1 -n 1 -v0 5054 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0

   2 species (with external immigration and death) example:  (OUTPUT_VARIABLES_GENUINE will be 5) 
   .~$ ./DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 2 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.1 -n 1 -v0 10105 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0
   Notice that the 2nd species (in green) does not die or externally immigrate.

   .~$ ./DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 2 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.5 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0
   Notice that the 2nd species (in green) does not die or externally immigrate.

   10 species (with external immigration and death) example (OUTPUT_VARIABLES_GENUINE will be 13) 
   .~$ ./DIFFUSION_S_RESOURCES -y0 1 -y2 1 -HS 10 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -H0 0.01 -H1 0.5 -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0

   Important notice: MODEL.h contains a 'define' of No_of_RESOURCES_MAXIMUM. If you overcome that
   limit, the program crashes. 

   4 species (4 different species ---R, A, RA, ARA---, therefore OUTPUT_VARIABLES_GENUINE is 7) 
   .~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100  -n 2 -v0 0 -v1 1 -G0 1 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1100.0

   MODEL = DIFFUSION_1R1C
   Single patch (-HM 1 -HX 1 -HY 1), and 3 species ---R, A, RA. Notice -H11 [Chi_C_0] -H12 [Eta_C_0]. If these two parameters are zero, no triplet formation, and the dynamics is equivalent to a 3D system, with only three local model variables.
   .~$ ./DIFFUSION_1R1C -y0 2 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 2 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 -tn 100 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -H1 0.0 -H6 0.5 -H9 10.0 -HK 2000  -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -H4 2.5 -H1 1.0 -H6 1.0 -H9 8.0 -H10 2.0 -H11 0.0 -H12 0.0

   MODEL = DIFFUSION_STOLLENBERG_3D
   Single patch (-HM 1 -HX 1 -HY 1), and 3 species ---R, A, RA. The dynamics do not include triplet formation. It is a 3D system, with only three local model variables (R, A, RA).
   .~$ ./DIFFUSION_STOLLENBERG_3D -y0 10 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 -tn 100 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -HK 2000 -H4 3.5 -H17 1.0 -H1 1.0 -H6 0.5 -H9 10.0 -H10 2.0

   .~$ ./DIFFUSION_STOLLENBERG_3D -y0 10 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 -tn 100 -t0 0.0 -t1 80.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 80.0 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 -HK 2000 -H4 1.5 -H17 1.0 -H1 1.0 -H6 0.5 -H9 10.0 -H10 10.0

   See denition_OutPut_Variables.c to understand the difference between Genuine Output Variables
   and plain model variables.

   MODEL = DIFFUSION_STOLLENBERG_4D
   Single patch (-HM 1 -HX 1 -HY 1), and 4 species ---RP, R, A, RA. The dynamics do not include triplet formation. It is a 4D system, with only four local model variables (RP, R, A, RA). 
   .~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 2000 -H20 20.0 -H1 1.0 -H3 5.0 -H6 0.5 -H9 20.0 -H10 2.0 -H4 5.0 -H17 1.0 
   This example produces a limit cycle. 
  
   .~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                                  -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 -tn 200 -t0 0.0 -t1 75.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 -G2 1 -G3 0.0 -G4 75.1 -G5 1 -G6 0.0 -G7 2000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 2000 -H20 20.0 -H1 1.0 -H3 5.0 -H6 0.5 -H9 5.0 -H10 2.0 -H4 5.0 -H17 1.0
   This example produces damped oscillations

   See denition_OutPut_Variables.c to understand the difference between Genuine Output Variables
   and plain model variables.
   
   -HuR -HuC are the jumping rates
   -H0 -H2 -H5  are the external immigration (Lambda_R_0, Lambda_R_1 and Lambda_C_0)
   -H20 is the establishment rate 
   -H1  -H3  -H6 are the death rates (Delta_R_0, Delta_R_1 for propagules, and Delta_C_0 for both searching and handling consumers)
   -H9  and -H10 are the Alpha_C_0 and Nu_C_0  Holling Type II model parameters 
   -H4  and -H17 are the production rates of propagules (Beta_R) and searching animals (Beta_C), respectively.  
   -xN  is INITIAL_TOTAL_POPULATION (involved in fixing Ini Conditions: include.Initial_Conditions.[].c)
   -HN  is No_of_INDIVIDUALS (involved in fixing some model parameters: see include.Parameter_Model.[].c )
   -xR  bool variable (0/1: Initial Conditions are not/yes re-scaled to the INITIAL_TOTAL_POPULATION value)

   Multi-patch system: a 100 x 100 grid (-HM 10000 -HX 100 -HY 100), and 4 species ---RP, R, A, RA. 
   The dynamics do not include triplet formation
   . ~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 \
                                   -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                                   -G2 1 -G3 0.0 -G4 75.1 -G5 1 -G6 0.0 -G7 2000000 \
                                   -tn 200 -t0 0.0 -t1 75.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 \
                                   -HuR 5.0 -HuC 1.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 200 -H20 20.0 \
                                   -H1 1.0 -H3 5.0 -H6 0.5 -H9 5.0 -H10 2.0 -H4 5.0 -H17 1.0

   . ~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 \
                                   -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                                   -tn 100 -t0 0.0 -t1 200.0 -t4 0 -tR 3 -xn 0 -xN 5.0 -HN 5.0 \
                                   -G2 1 -G3 0.0 -G4 200.1 -G5 1 -G6 0.0 -G7 200000 \
                                   -HuR 5.0 -HuC 1.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 20 -H20 20.0 \
                                   -H1 1.0 -H3 5.0 -H6 0.5 -H9 20.0 -H10 2.0 -H4 5.0 -H17 1.0

    . ~$ ./STOLLENBERG_4D_3 -y0 15 -y2 1 -HS 1 -HM 1600 -HX 40 -HY 40 \
                            -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                            -tn 100 -t0 0.0 -t1 200.0 -t4 0 -tR 3 -xn 0 -xN 5.0 -HN 5.0 \
                            -G2 1 -G3 0.0 -G4 200.1 -G5 1 -G6 0.0 -G7 32000 \
                            -HuR 5.0 -HuC 1.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 20 -H20 20.0 \
                            -H1 1.0 -H3 5.0 -H6 0.5 -H9 10.0 -H10 2.0 -H4 5.0 -H17 1.0

    . $ ./STOLLENBERG_4D_3  -y0 15 -y2 1 -HS 1 -HM 6400 -HX 80 -HY 80 \
                            -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2  \
                            -tn 400 -t0 0.0 -t1 600.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                            -G2 1 -G3 0.0 -G4 600.1 -G5 1 -G6 0.0 -G7 128000 \
                            -HuR 5.0 -HuC 1.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 20 -H20 20.0 \
                            -H1 1.0 -H3 5.0 -H6 0.5 -H9 15.5 -H10 2.0 -H4 5.0 -H17 1.0

   . ~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                                   -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \ 
                                   -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 \ 
                                   -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 2000 \
                                   -HuR 0.0 -HuC 0.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 2000 -H20 20.0 \
                                   -H1 1.0 -H3 5.0 -H6 0.5 -H9 20.0 -H10 2.0 -H4 5.0 -H17 1.0 
   
   . ~$ ./DIFFUSION_STOLLENBERG_4D -y0 15 -y2 1 -HS 1 -HM 100 -HX 10 -HY 10 \
                                   -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                                   -G2 1 -G3 0.0 -G4 200.1 -G5 1 -G6 0.0 -G7 2000 \
                                   -tn 200 -t0 0.0 -t1 200.0 -t4 0 -tR 10 -xn 0 -xN 5.0 -HN 5.0 \
                                   -HuR 5.0 -HuC 1.0 -H0 0.0 -H2 0.0 -H5 0.0 -HK 20 -H20 20.0 \
                                   -H1 1.0 -H3 5.0 -H6 0.5 -H9 5.0 -H10 2.0 -H4 5.0 -H17 1.0

  . ~$ ./DIFFUSION_AZTECA_4D -y0 17 -y2 1 -HS 1 -HM 6400 -HX 80 -HY 80 \
                             -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2  \
                             -tn 400 -t0 0.0 -t1 600.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                             -G2 1 -G3 0.0 -G4 600.1 -G5 1 -G6 0.0 -G7 128000 \
                             -HuR 1.0 -HuC 5.0 -H0 0.0 -H5 0.0 \
                             -HK 100 -H7 5 \
                             -H1 5.0 -H3 0.5 -H6 2.5 -H8 10.0 \
                             -H9 15.5 -H10 2.0 -H4 5.0 -H20 20.0

  ~$ ./DIFFUSION_AZTECA_4D -y0 17 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                           -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                           -tn 100 -t0 0.0 -t1 100.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 \
                           -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 3000 \
                           -HuR 5.0 -HuC 1.0 -H0 0.0 -H5 0.0 \
                           -HK 1000 -H7 3 \
                           -H1 5.0 -H3 0.5 -H6 2.5 -H8 10.0 \
                           -H9 15.5 -H10 2.0 -H4 5.0 -H20 20.0

  ~$ ./DIFFUSION_AZTECA_4D -y0 17 -y2 1 -HS 1 -HM 6400 -HX 80 -HY 80 \
                           -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2  \
                           -tn 400 -t0 0.0 -t1 600.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                           -G2 1 -G3 0.0 -G4 600.1 -G5 1 -G6 0.0 -G7 128000 \
                           -HuR 1.0 -HuC 5.0 -H0 0.0001 -H5 0.0001 \
                           -HK 1000 -H7 100 \
                           -H1 5.0 -H3 0.5 -H6 2.5 -H8 10.0 \
                           -H9 10000.5 -H10 2.0 -H4 10.0 -H20 10.0 -tE 2.1

  For the AZTECA_4D model,
  -HuR and -HuC are jumping rates of workers and flies, respectively.
  -H0 and -H5 are external immigrations rates (Lambda_R_0 and Lambda_C_0). Usually equal to zero
  -H7 is the maximum number of nests (queens) per local area, cell or local community (Lambda_C_1 !!!)
  -HK is the maximum number of workers per nest. 
  -H20 is the establishment rate (of a new nest, Eta_R)
  -H4 is the per capita rate at which queens produce new workers (Beta_R)
  -H9 is the attack rate (Alpha_C_0)
  -H10 is the larval development rate (Nu_C_0)  
  -H1 -H3 -H6 -H8 are death rates (workers, queens, flies, and parasitezed workers, respectively) 
  -xN n_0. For instance, the initial condition is (0, n_0, 0, 0) at the center of the grid and 
           (0, n_0, n_0, n_0) at a number of randomly chosen cells across the grid. 

  AZTECA_4D_0 only uses one single carrying capacity, the total number of pontential nesting trees, 
  which can be initialized by using -HK [Number] as an input argument 
  ~$ ./DIFFUSION_AZTECA_4D_0 -y0 18 -y2 1 -HS 1 -HM 6400 -HX 80 -HY 80 \
                             -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2  \
                             -tn 100 -t0 0.0 -t1 60.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                             -G2 1 -G3 0.0 -G4 60.1 -G5 1 -G6 0.0 -G7 128000 \
                             -HuR 1.0 -HuC 5.0 -H0 0.0001 -H5 0.0001 \
                             -HK 10 \
                             -H1 5.0 -H3 0.5 -H6 2.5 -H8 10.0 \
                             -H9 100.5 -H10 2.0 -H4 10.0 -H20 10.0 -tE 2.1

  ~$ ./DIFFUSION_AZTECA_4D_0 -y0 18 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                            -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                            -tn 100 -t0 0.0 -t1 100.0 -t4 0 -tR 10 -xn 0 -xN 500.0 -HN 500.0 \
                            -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 3000 \
                            -HuR 5.0 -HuC 1.0 -H0 0.0 -H5 0.0 \
                            -HK 1000 \
                            -H1 5.0 -H3 0.5 -H6 2.5 -H8 10.0 \
                            -H9 15.5 -H10 2.0 -H4 5.0 -H20 20.0

  For the AZTECA_4D_0 model,
  -HuR and -HuC are jumping rates of workers and flies, respectively.
  -H0 and -H5 are external immigrations rates (Lambda_R_0 and Lambda_C_0). Usually equal to zero
  -HK is the maximum number of potential nesting trees (queens) per local area, cell or local community (K_R)  
  -H20 is the establishment rate (of a new nest, Eta_R)
  -H4 is the per capita rate at which queens produce new workers (Beta_R)
  -H9 is the attack rate (Alpha_C_0)
  -H10 is the larval development rate (Nu_C_0)  
  -H1 -H3 -H6 -H8 are death rates (workers, queens, flies, and parasitezed workers, respectively) 
  -xN n_0. For instance, the initial condition is (0, n_0, 0, 0) at the center of the grid and 
           (0, n_0, n_0, n_0) at a number of randomly chosen cells across the grid. 

  ~$ ./DIFFUSION_AZTECA_4D_1 -y0 19 -y2 1 -HS 1 -HM 6400 -HX 80 -HY 80  \
                        -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2 \
                        -tn 300 -t0 0.0 -t1 300.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                        -G2 1 -G3 0.0 -G4 300.0 -G5 1 -G6 0.0 -G7 50000 \
                        -HuR 1.0 -HuQ 2.5 -HuC 5.0 \
                        -H0 0.0001 -H2 0.0001 -H5 0.0001 \
                        -HK 10.0 \
                        -H1 1.5 -H3 1.0 -H6 1.5 -H8 1.0 \
                        -H9 50.0 -H10 2.0 -H4 8.5 -H20 1.0

  For the AZTECA_4D_1 model,
  -HuR -HuQ and -HuC are jumping rates of workers, queens and flies, respectively.
  -H0 -H2 and -H5 are external immigrations rates (Lambda_R_0, Lambda_R_1, and Lambda_C_0). Usually equal to zero
  -HK is the maximum number of potential nesting queens per local area, cell or local community (K_R)  
  -H20 is the food per unit worker carried back to the nest (times a conversion factor) (Eta_R)
  -H4 is the per capita rate at which queens produce new workers (Beta_R)
  -H9 is the attack rate (Alpha_C_0)
  -H10 is the larval development rate (Nu_C_0)  
  -H1 -H3 -H6 -H8 are death rates (workers, queens, flies, and parasitezed workers, respectively) 
  -xN n_0. For instance, the initial condition is (0, n_0, 0, 0) at the center of the grid and 
           (0, n_0, n_0, n_0) at a number of randomly chosen cells across the grid.                       

  Tradeoff: Between Eta and Beta based on an invasion criterion (r_0): 
  Table->Beta_AP[i] = r_0 * (Table->Eta_RP[i] + Table->Delta_RP[i])/(Table->Eta_RP[i]*(1.0-2.0*p)) * Table->Delta_AP[i];

  ~$ ./DIFFUSION_ECOEVO_PLANTS  -y0 20 -y2 1 -HS 100 -HM 1 -HX 1 -HY 1  \
                                -n 1 -G0 1 -G1 1 \
                                -tn 300 -t0 0.0 -t1 300.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                                -G2 1 -G3 0.0 -G4 300.0 -G5 1 -G6 0.0 -G7 100 \
                                -HuR 0.0 \
                                -H0 0.0 \
                                -HK 1000.0 \
                                -H1 1.5 -H3 1.0 -H6 0.25 -H8 4.0 \
                                -Hp1 0.0001 \
                                -Hp2 2.0 \
                                -H20 5.0 -H4 2.0
  
  Tradeoff: Between Eta and Beta based on the amount of free space (z_0): 
  Table->Beta_AP[i] = z_0*(Table->Delta_RP[i]+Table->Eta_RP[i])/( z_0*(Table->Delta_RP[i]+Table->Eta_RP[i])-Table->Delta_RP[i]) * Table->Delta_AP[i]; 
 
  ~$ ./DIFFUSION_ECOEVO_PLANTS  -y0 20 -y2 1 -HS 100 -HM 1 -HX 1 -HY 1  \
                                -n 1 -G0 1 -G1 1 \
                                -tn 300 -t0 0.0 -t1 300.0 -t4 0 -tR 1 -xn 0 -xN 5.0 -HN 5.0 \
                                -G2 1 -G3 0.0 -G4 300.0 -G5 1 -G6 0.0 -G7 1000 \
                                -HuR 0.0 \
                                -H0 0.0 \
                                -HK 1000.0 \
                                -H1 1.0 -H3 0.5 -H6 0.5 -H8 8.0 \
                                -Hp1 0.0001 \
                                -Hp2 0.8 \
                                -H20 5.0 -H4 2.0

  In case you want to represent local variables across cells, keep in mind that there are 2 local variables (RP and R)
  per species. If -HS 100, then you have 200 local variables to represent. For example, you can decide to plot the
  following 16 sp in a 4 x 4 subpanel organization:

  -G0 4 -G1 4 \
  -n 16 -v0 69 -v1 71 -v2 73 -v3 75 -v4 77 -v5 79 -v6 81 -v7 83 -v8 85 -v9 87 -v10 89 -v11 91 -v12 93 -v13 95 -v14 97 -v15 99

  Other subplotting: -n 4 -v0 0 -v1 1 -v2 2 -v3 3 -G0 2 -G1 2. Subplotting is called from int M_O_D_E_L( Parameter_Table * Table )
  in MODEL.c. 
  
  Otherwise:
          -n 1, and only one plot controls de the deterministic evolution of the different types in a single bar plot
                  
  -HN (No_of_INDIVIDUALS) is not in use in ECOEVO_PLANTS model (only -xN, INITIAL_TOTAL_POPULATION )
                       
  If -HM 1, there is only one patch, then -HuR (or -Hu) should be 0.0. Propagule movement between patches does not make sense.

  For the ECOEVO_PLANTS model,
  -HS is the number of resource type (discretization of the trait axis) (No_of_RESOURCES)
  -HuR (or -Hu) is the jumping rate of resource propagules.
  -H0 is the external immigrations rate (Lambda_R_0) of resource propagules. Usually equal to zero
  -HK is the maximum number of potential sites per local area, cell or local community (K_R)  
  -H20 is a establishment rate of resource propagules into adult plants (Eta_R)
  -H4 is the per capita rate at which adult plants produce new propagules (Beta_R)
  -H1 -H3 are death rates of resource probagules (Delta_R_0) and adult plants (Delta_R_1), respectively) 
  -H6 -H8 are max and min parameter values defining a tradeoff (Delta_C_0 is the min and Delta_C_1 is the MAX)
  -Hp1 is the mutation rate (p_1)    0 < p_1 < 1
  -Hp2 is the tradeoff factor (p_2)  0 < p_2  (no greater than 10, because there is an assert(p_2 < 10.0) )    
  -xN n_0. For instance, the initial condition is (n_0, 0) for species 0 at the center of the grid and 
           (n_0, 0) at a number of randomly chosen cells across the grid.

   MacArthur and Rosenzweig (two species 3D, R, A, RA):
   .~$ ./DIFFUSION_MR -y0 7 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                      -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \
                      -tn 100 -t0 0.0 -t1 40.0 -t4 0 -tR 10 -xn 0 -xN 50.0 -HN 50.0 \ 
                      -G2 1 -G3 0.0 -G4 40.0 -G5 1 -G6 0.0 -G7 2000 \
                      -H1 0.0 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 \
                      -H4 25.0 -H1 0.0 -H6 5.0 -H9 17.0 -H10 5.0 -HK 2000.0

   Coarse-grained 2D system (two species 1R and 1C):
   .~$ ./DIFFUSION_1R1C_2D -y0 4 -y2 1 -HS 1 -HM 4 -HX 2 -HY 2 \ 
                           -n 2 -v0 0 -v1 1 -G0 1 -G1 2 \
                           -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 4 -xn 0 -xN 50 -HN 50 \ 
                           -G2 1 -G3 0.0 -G4 40.0 -G5 1 -G6 0.0 -G7 100.0 \
                           -H0 0.01 -H5 0.01 -H6 0.5

   MODEL = DIFFUSION_BD_2D
   Feeding experiments at a constant number of total consumers: 
   .~$ ./DIFFUSION_BD_2D -y0 13 -y2 1 -HS 1 -HM 1 -HX 1 -HY 1 \
                         -n 2 -v0 0 -v1 1 -G0 1 -G1 2 \
                         -tn 50 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 \
                         -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 20 \
                         -HK 10000 -HuR 0.0 -HuC 0.0 -H0 0.0 -H5 0.0 \
                         -H9 2.5 -H10 10.0 -H11 100.0 -H12 1.0 -Hp1 0.3725 -Hp2 0.5 -HN 20

   MODEL = DIFFUSION_HII_nD
   Feeding experiments at a constant number of total consumers (-HN 20): 
   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \ 
                          -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \ 
                          -tn 50 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 \
                          -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 20 \
                          -HK 10000 -HuR 0.0 -HuC 0.0 -H0 5.0 -H2 1.0 -H5 0.0 \
                          -H9 2.5 -H10 10.0 -Hp1 0.3725 -Hp2 0.5 -HN 20

   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \
                                 -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \
                                 -tn 10 -t0 0.0 -t1 2.5 -t4 0 -tR 100 -xn 0 -xN 20.0 \
                                 -G2 1 -G3 0.0 -G4 2.5 -G5 1 -G6 0.0 -G7 8 \
                                 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 5.0 -H2 1.0 -H5 0.0 \
                                 -H9 10.5 -H10 10.0 -Hp1 0.3725 -Hp2 1.0 -HN 20 -tE 0.2

   .~S ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \
                          -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \
                          -tn 10 -t0 0.0 -t1 2.5 -t4 0 -tR 100 -xn 0 -xN 20.0 \
                          -G2 1 -G3 0.0 -G4 2.5 -G5 1 -G6 0.0 -G7 20 \
                          -HK 10000 -HuR 0.0 -HuC 0.0 -H0 5.0 -H2 1.0 -H5 0.0 \
                          -H9 10.5 -H10 0.1  -Hp1 0.3725 -Hp2 1.0 -HN 20 -tE 0.2
   
Relevant input arguments for model DIFFUSION_HII_nD:
   -HK  [N or Total Carrying Capacity (for resources)]
   -Hp1 [approx y_R = f_i * p1 * N: Different resource levels at time 0.0, with f_i = 0.5, 0.3, 0.2, but this can be changed] 
   -Hp2 [Total No of Free Consumers at time 0.0: p2 * No_of_CONSUMERS]
   -HN 20 [Total No of CONSUMERS]
   -H9 and -H10 [Alpha_C_0 and Nu_C_0 for 1st Resource Type. Alpha_C_0 is the same for all resource types, but this can be changed.  
   -H0 and -H2 are Lambda_R_0 and Lambda_R1, which are overloaded to create extra handling rates (Nu), for the 2nd and 3rd resource type. 
   Notice that for this model -H5 should be always zero (No immigration of consumers: No of CONSUMERS is constant). 
   -H11 and H12 [Xhi_C_0 and Eta_C_0 are not relevant here. They are only in models with predator interference].
   
   MODEL = DIFFUSION_ECO_1B1P (-y2 22)  
  .~$ ./DIFFUSION_ECO_1B1P -y0 22 -y2 1 -HS 1 -HN 1 -HM 1 -HX 1 -HY 1 \
                              -n 2 -v0 0 -v1 1 -G0 1 -G1 2 \
                              -G2 1 -G3 0.0 -G4 10.0 -G5 1 -G6 0.0 -G7 1000 \
                              -tn 100 -t0 0.0 -t1 10.0 -t4 0 -tR 10 -xn 0 -xN 50.0  \
                              -HK 1000 -H4 2.5 -H1 1.0 -H3 1.0 -H2 1.0\
                              -H9 0.5 -H10 0.5 \
                              -Hp1 0.01 -Hp2 0.5 -H11 0.3 \
                              -HuR 0.0 -HuC 0.0 -H0 0.0

  .~$ ./ECO_1B1P_4 -y0 22 -y2 1 -HS 1 -HN 1 -HM 1 -HX 1 -HY 1 \
                   -n 2 -v0 0 -v1 1 -G0 1 -G1 2  \
                   -G2 1 -G3 0.0 -G4 100.0 -G5 1 -G6 0.0 -G7 200 \
                   -tn 100 -t0 0.0 -t1 100.0 -t4 0 -tR 10 -xn 0 -xN 50.0 \
                   -HK 500 -H4 5.5 -H1 0.1 -H3 0.1 -H2 10.0 \
                   -H9 0.1 -H10 0.2 \
                   -Hp1 0.01 -Hp2 0.005 -H11 0.00025 \
                   -HuR 0.0 -HuC 0.0 -H0 0.0
                              
  .~$ ./ECO_1B1P_4 -y0 22 -y2 1 -HS 1 -HN 1 -HM 2500 -HX 50 -HY 50 \
                     -n 2 -v0 0 -v1 1 -G0 1 -G1 2 \
                     -G2 1 -G3 0.0 -G4 1000.0 -G5 1 -G6 0.0 -G7 300 \
                     -tn 500 -t0 0.0 -t1 1000.0 -t4 0 -tR 20 -xn 0 -xN 40.0 \
                     -HK 100 -H4 2.5 -H1 1.0 -H3 1.0 -H2 1.0 \
                     -H9 0.25 -H10 0.7 \
                     -Hp1 0.01 -Hp2 0.17 -H11 0.51 \
                     -HuR 2.0 -HuC 4.0 -H0 0.0

  (see COMPILATION_OPTIMIZATION_LEVELS.sh for compilation optiopns)               
  Relevant input arguments for model DIFFUSION_1B1P:  
   -HN 1 [No_of_PLASMIDS]: always 1 plasmid (for this model, 1B1P) (As input arg, coded as No_of_INDIVIDUALS) 
   -HS 1 [No_of_STRAINS]:  always 1 strain (for this model, 1B1P) (As input arg, coded as No_of_RESOURCES)
   -HM 1 [No_of_PATCHES]:  it could be more than one, but keep it between 1 and 10000 maximum
   -HX 1 [No_of_CELLS_X]:  it could be more than one, but keep it between 1 and 100 maximum, accordingly to the number of patches
   -HY 1 [No_of_CELLS_Y]:  it could be more than one, but keep it between 1 and 100 maximum, accordingly to the number of patches
   -HK  Total Carrying Capacity in every local population
   -H1:  Delta_R_0,  basal per capita death rate for both plasmid-free and plasmid-carrying cells                                                                               
   -H3:  Delta_R_1,  stress-induced per capita death rate for both plasmid-free and plasmid-carrying cells
   -H4:  Beta_R, cell reproduction rate
   -Hp1 Segregation error at reproduction
   -Hp2 Competition-induced Death Probability
   -H11: Plasmid Transmission Probability  (coded as Chi_C_0, but labeled as $\Xhi$ in the paper)
   -H9 Alpha_C_0, Reproduction Cost (associated to the presence of the plasmid) 
   -H10 Nu_C_0, Resistance to stress-induced death (associated to the presence of a plasmid)    
   -H0 Lambda_R_0 is the external immigration rate for both plasmid free and plasmid carrying cells.
   -H2 Lambda_R_1  Conjugation/encounter rate (coded as Lambda_R_1, but labeled as $\gamma_0$ in the paper)
   -HuR -HuC are the jumping rates for plasmid free and plasmid carrying cells (only relevant if there are more than one cell or patch in the system)
   -y2 Type of network: 1 for a 2D lattice, with periodic boundary conditions.
   -H12 Eta_C_0 is not relevant here. They is only in models with predator interference.
   -H5 is not relevant here. (Lambda_C_0, an external immigration of a second species,
       but both free plamsid and carrying cells have the same external immigration rate, Lambda_R_0)
   */

int main(int argc, char **argv)
{
  int i;
  Parameter_Table Table;
  Time_Control Time;
  Time_Dependence_Control Time_Dependence;
  P_ARG = &Table;

#include "default.c"

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

#include <include.Output_Variables.default.aux.c>
  
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C(   &Table );
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( &Table, Index_Output_Variables );
  printf(" Parameter_Table structure has been correctly allocated and initiated\n");

  /* B E G I N : Reserving memmory for Parameter Space */
#include <include.Parameter_Space.default.aux.c>
  if( No_of_PARAMETERS == Table.TOTAL_No_of_MODEL_PARAMETERS ) {
    /* Full parameter space is in place. See also Model_Variables_Code.c */
    for(i=0; i<Table.TOTAL_No_of_MODEL_PARAMETERS; i++) Index[i] = Table.Index[i];
    No_of_PARAMETERS = Table.TOTAL_No_of_MODEL_PARAMETERS;
  }
  Parameter_Space * Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  Parameter_Space_Alloc( Space, No_of_PARAMETERS, d);
  Parameter_Space_Initialization( Space, No_of_PARAMETERS, TOLERANCE, MAX_No_of_ITERATIONS,
    d, Index, Ranges, Acc);
  Table.S = Space;
  printf(" Parameter_Space structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */

#include <gsl_random_number_Setup.c>
  // #if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); Print_Press_Key(1,0,".");
  /*   END: Checking Random Number Generator Setup */
  // #endif

  if (TYPE_of_TIME_DEPENDENCE == 0) {
    printf(" Time_Control structure will be allocated: \n");
    printf(" %d output variables of length %d points will be allocated\n",
	   SUB_OUTPUT_VARIABLES, I_Time);
    T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( &Time, &Table, I_Time);
    T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D( &Time, &Table, I_Time);
    printf(" Time_Control structure has been correctly allocated and set up\n");
  }
  else {
    #include <include.Time_Dependence_Control.default.aux.c>
    printf(" Time_Dependence_Control and Time_Control structures will be allocated: \n");
    printf(" %d output variables of length %d points will be allocated\n",
    SUB_OUTPUT_VARIABLES, I_Time);
    Time_Dependence_Control_Alloc(&Time, &Time_Dependence, &Table,
				  I_Time, TIME_DEPENDENT_PARAMETERS, No_of_COVARIATES);

    int No_of_EMPIRICAL_TIMES = I_Time;
    // Number of columns in the data files of time-dependent parameters
    Time_Dependence_Control_Upload(&Time, &Time_Dependence, &Table,
				   I_Time, No_of_EMPIRICAL_TIMES,
				   TIME_DEPENDENT_PARAMETERS, TYPE_of_TIME_DEPENDENCE,
				   TYPE_0_PARAMETERS, TYPE_1_PARAMETERS, TYPE_2_PARAMETERS,
				   No_of_COVARIATES,
				   dependent_parameter, forcing_pattern,
				   "File_of_Covariates.dat", Name_of_FILE[0] );
    printf(" Both Time_Control and Time_Dependence_Control structures have been\n");
    printf(" correctly allocated and set up\n");
  }

#if defined CPGPLOT_REPRESENTATION
  Table.CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  //Table.CPG_STO = A_C_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T (1, SUB_OUTPUT_VARIABLES, I_Time, 0, "/TPNG");
  Table.CPG_STO = A_C_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T (1, SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  
  printf(" Two Parameterh_CPGPLOT plotting structures have been correctly allocated and initiated\n");
  printf(" These will open two windows (or two ploting devices of the same kind)\n");
  printf(" Table.CPG will store deterministic dynamic variables to plot\n");
  printf(" Table.CPG_STO will store stochastic dynamic variables to plot\n");
  printf(" As a consquence, deterministic and stochastic dynamics can be plotted\n");
  printf(" on the same device to compare (as it is done here, indicated by the first\n");
  printf(" input argument (0) of the A_Ch_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T function).\n");
  printf(" Alternatively, two different devices (two different pdf files, for instance)\n");
  printf(" can be used, if required (1).\n");
#endif

  if(Table.TYPE_of_MODEL == 12 || Table.TYPE_of_MODEL == 13 || Table.TYPE_of_MODEL == 14) 
    /* Models where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);

    /* Models where the TOTAL_No_of_CONSUMERS is not a CONSTANT 4D Consumer-Resource Models */
    /* MODEL = DIFFUSION_AZTECA_4D  and MODEL = DIFFUSION_STOLLENBERG_4D                    */
    // if(Table.TYPE_of_MODEL == 15 || Table.TYPE_of_MODEL == 17) 
    // Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table); //Under construction

  if(Table.TYPE_of_MODEL == 16) { // DIFFUSION_HII_nD
    /* Also model where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    /* and they feed on multiple resources                      */
    Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
    Resetting_Alpha_Nu_Vectors (&Table);
    Resetting_Multiresource_Levels (&Table);  
    Writing_Alpha_Nu_Theta_Vectors(&Table);  
  }

#ifndef DIFFUSION_ECO_PLASMIDS /* First Test without Deterministic Dynamics */     
    Parameter_Values_into_Parameter_Table(&Table);

  /* Deterministic Time Dynamics */  
  if(Table.TYPE_of_MODEL == 20) // DIFFUSION_ECOEVO_PLANTS
    Ressetting_Species_Characteristic_Parameters (&Table);

  M_O_D_E_L( &Table );
#endif

  // Some models (such as DIFFUSION_1R1C_2D and so on) do not have a stochastic
  // counter-part implemented yet!

  if (Table.TYPE_of_MODEL == 21 || Table.TYPE_of_MODEL == 17 || Table.TYPE_of_MODEL == 18 || 
      Table.TYPE_of_MODEL == 19 || Table.TYPE_of_MODEL == 20 || Table.TYPE_of_MODEL == 2  || 
      Table.TYPE_of_MODEL == 8  || Table.TYPE_of_MODEL == 10 || Table.TYPE_of_MODEL == 15 || 
      Table.TYPE_of_MODEL == 12 || Table.TYPE_of_MODEL == 13 || Table.TYPE_of_MODEL == 14 || 
      Table.TYPE_of_MODEL == 9  || Table.TYPE_of_MODEL == 16 || Table.TYPE_of_MODEL == 22){
    
    M_O_D_E_L___S_T_O( &Table );
  }
  else{
    Print_Press_Key(1,1,"Stochastic optimization has not been implemented for this model.\n");
  }
  
  /* BEGIN : -------------------------------------------------------------------------
   */
  char boundary_File[80];
  sprintf(boundary_File, "boundary_Model_Parameter.c");
  write_Parameter_Table___RANGES___VALUES___LATEX ( "Latex_Parameter_Table.tex",
                                                    boundary_File,
                                                    &Table,
                                                    Space->P_MAX->data,
                                                    Space->P_min->data, Space->No_of_PARAMETERS );
  /*  END : ------------------------------------------------------------------------*/

  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */
#if defined CPGPLOT_REPRESENTATION
  #include <include.CPG.default.free.c>
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG_STO, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif

#include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );

#include <include.Initial_Conditions.default.free.c>

#include <include.Output_Variables.default.free.c>

#include <include.Time_Dependence_Control.default.free.c>
  if (TYPE_of_TIME_DEPENDENCE == 0) T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( &Time, &Table );
  else                        Time_Dependence_Control_Free( &Time_Dependence, &Table );

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */

  printf("\nEnd of progam\n");
  return (0);
}
  
