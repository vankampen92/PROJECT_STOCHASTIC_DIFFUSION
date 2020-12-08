#include "GSL_stat.h"

/* Discrete Sampling provide a set of functions required in 
   a lot of stochastic models, one step stochastic simulations and 
   so on. 
   These functions are a GSL alternative to the functions in
   stat_Discrete_Sampling.c. Obviously, all functions are 
   overloaded and the compilation of both set of functions
   at the same time is forbidden. It would produce an ERROR. 

   Total independence from numerical recipes !!! 
*/

extern gsl_rng * r; /* Global generator defined in main.c */

#define RANDOM gsl_rng_uniform_pos(r)

int Discrete_Sampling_Rejection_Method(double p_max, double * a, int NoEvents)
{
  int stat_Bool;
  int xr;
  double R_A_T_E;

  stat_Bool = 0;
  while(stat_Bool == 0){

    xr = 1 + gsl_rng_uniform_int(r, NoEvents);
    // xr = (int)random_int_lrand48(1, NoEvents);
    R_A_T_E = a[xr-1]/p_max;
    if(RANDOM < R_A_T_E)  stat_Bool = 1;
  }  
  
  return xr;
}

int Discrete_Sampling_Old(double * a, int NoEvents)
{
  int j, kount;
  double xr, rate;
  double * R_A_T_E;
  
  xr = RANDOM;
  R_A_T_E = (double *)calloc(NoEvents+1, sizeof(double) );
  R_A_T_E[0] = 0.; 
  for (j=0; j<NoEvents; j++)
    R_A_T_E[j+1]= R_A_T_E[j] + a[j];
  rate = R_A_T_E[NoEvents];
  for (j=0; j<NoEvents; j++) 
    R_A_T_E[j+1]/=rate; /* Normalizing rates...*/ 
  kount = 0;
  for (j=0; j<NoEvents; j++)               /* Deciding elemetary event to occur */
    if((xr>R_A_T_E[j]) && (xr<=R_A_T_E[j+1])) kount = j+1;
  
  free (R_A_T_E);
  return kount;
}

int Discrete_Sampling(double *a, int NoEvents)
{
  int j, kount;
  double * R_A_T_E= (double *)calloc( NoEvents, sizeof(double) );

  R_A_T_E[0] = a[0];
  for (j=1; j<NoEvents; j++)
    R_A_T_E[j]= R_A_T_E[j-1] + a[j];
  
  kount = Discret_Sampling_High_Performance(R_A_T_E[NoEvents-1], R_A_T_E, NoEvents);
 
  free(R_A_T_E);
  return kount;
}

int Discrete_Sampling_Cummulative(double rate, double *R_A_T_E, int NoEvents)
{
  /* . rate is the normalization constant.
     . r[] must store the NON-NORMALIZED cummulative distribution probability,
     where $rate$ is the normalizing factor. If rate=1., the algorithm could be
     improved but actually works quite well too.
     . If r[] stores the discret probabililty distribution function, this algorithm 
     does not work at all. Then, the Discret_Sampling() algorith must be used instead.
  */
  int j, kount;
  double x_1, x_2, xr;
  
  xr = RANDOM;
  kount = 1;
  x_2 = 0.;
  for (j=0; j<NoEvents; j++){/* Deciding elemetary event to occur */
    x_1 = x_2;
    x_2 = R_A_T_E[j]/rate;
    if((xr > x_1) && (xr <= x_2)) {
      kount = j+1;
      return kount;
    }
  }
  return(0);
}

int Discret_Sampling_High_Performance(double rate, double *R_A_T_E, int NoEvents)
{
  /* . rate is the normalization constant.
     . R_A_T_E[] must store the non-normalized cummulative probability distribution,
     where rate is the normalizing factor. If rate=1., the algorithm could be
     improved but actually works quite well too.
     . If R_A_T_E[] stores the discret probabililty distribution function, 
     rather than the cummulative probability distribution, this algorithm 
     does not work at all. Then, the Discret_Sampling() algorith must be used instead.
  */
  int j1,j2,jm, kount;
  int stat_Bool;
  double x_1, x_2, x_m, xr;

  if(NoEvents == 1) 
    return(1);
 
  xr = RANDOM;

  if(NoEvents == 2){
    x_1 = 0.; x_2 = R_A_T_E[0]/rate;
    if((xr > x_1) && (xr <= x_2))
      return(1);
    else 
      return(2);
  }
  else if(NoEvents > 2 && NoEvents < 10){
    kount = Discrete_Sampling_Cummulative(rate, R_A_T_E, NoEvents);
    return (kount);
  }
  else{
    x_1 = 0.; x_2 = R_A_T_E[0]/rate;
    if((xr > x_1) && (xr <= x_2))
      return(1);
    x_1 = R_A_T_E[NoEvents-2]/rate; x_2 = R_A_T_E[NoEvents-1]/rate;
    if((xr > x_1) && (xr <= x_2))
      return(NoEvents);

    stat_Bool = 0;
    j1 = 0;            x_1 = R_A_T_E[j1]/rate;
    j2 = NoEvents-1;   x_2 = R_A_T_E[j2]/rate;
    while(stat_Bool == 0){ 
      jm = (j1 + j2)/2; /* Integer division by 2 */
      x_m = R_A_T_E[jm]/rate;
      if(j1 == j2-1)  
	stat_Bool = 1;       /* Brackening has succedeed  p(j1)< xr <= p(j2), p(i) = R_A_T_E[i]/rate */
      else{
	if((xr > x_m) && (xr <= x_2)){
	  x_1 = x_m;
	  j1  = jm;
	}
	else{
	  x_2 = x_m;
	  j2  = jm;
	}
      }
    }
  }
  return (j2+1);
}

