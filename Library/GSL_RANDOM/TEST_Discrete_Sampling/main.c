#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Discrete_Sampling.h"

int main()
{
    int i; 
    int M; 
    int Event; 

    M  = 100000; /* M plays the role of the number of events. */

    printf("Creating, first, a discrete distribution of M values:\n");
    printf("Enter and integer-valued number, M, please... ");
    scanf("%d", &M); 

    double * Rates = (double *)calloc(M, sizeof(double));
    
    for (i=0; i<M; i++) {
        Rates[i] = drand48();
        printf( "Rate[%d] = %g\n", i, Rates[i] );
    }    

    Discrete_Sampling(Rates, &M, &Event);

    printf(" The winner event is... the %d-th!!!\n", Event);
    printf(" (from 1 to %d)\n", M);

    free(Rates);
    return(0);
}