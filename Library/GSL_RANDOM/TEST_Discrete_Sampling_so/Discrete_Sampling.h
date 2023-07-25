/* Discrete Sampling provide a set of functions required in 
   a lot of stochastic models, one step stochastic simulations and 
   so on. 
*/
int Discrete_Sampling(double *a, int NoEvents);
int Discret_Sampling_High_Performance(double rate, double *R_A_T_E, int NoEvents);
int Discrete_Sampling_Cummulative(double rate, double *R_A_T_E, int NoEvents);