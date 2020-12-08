#include "GSL_stat.h"

void GSL_Init_Random_Seed(const gsl_rng * r)
{
  /* This function seed the GSL Random Generator r
     with a seed which is different for each initialization
     according to current computer time 
  */
        unsigned long int     seed;
        time_t  nowtime;
        struct  tm *preztime;

        time(&nowtime);
        preztime = localtime(&nowtime);
        seed = (int)((preztime->tm_sec+1)*(preztime->tm_min+1)*
                (preztime->tm_hour+1)*(preztime->tm_year)*(preztime->tm_year));
        if(seed%2==0) seed++;

	printf(" Random Number Seed: %lu\n", seed);

        gsl_rng_set(r, seed);
}

void GSL_Init_Random_Seed_from_File(const gsl_rng * r)
{
        unsigned long int     seed;

	//seed = 100;       
	/* A script have changed the environmetal variable
	   GSL_RNG_SEED before the execution of the program
	   starts. This value has been set to gsl_rng_default_seed
	   in gsl_random_number_Setup.c. This setup is 
	   always done at the start of code execution.
	*/

  	printf ("GSL_RNG_SEED = %lu\n", gsl_rng_default_seed);

	gsl_rng_set(r, gsl_rng_default_seed);

	printf ("first value = %lu\n", gsl_rng_get (r));
}

