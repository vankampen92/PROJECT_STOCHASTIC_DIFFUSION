#include <MODEL.h>

extern gsl_rng * r; /* Global generator defined in main.c */

void show_DoubleMatrix(double **M, int Nx, int Ny)
{
  int i,j, n_sys;

  /*n_sys = system("clear");*/
  printf("Printing matrix...\n");

  printf("* :\t");
  for(i=0; i<Ny; i++) printf("%4d  ",i);
  printf("\n");
  printf("   \t");
  for(i=0; i<Ny; i++) printf("..... ");
  printf("\n");

  for(i=0; i<Nx; i++){
    printf("%d :\t",i);
    for(j=0; j<Ny; j++) printf("%4.3g  ",M[i][j]);
    printf("\n");
  }
  printf("\n\n");
  /* getchar(); */
}