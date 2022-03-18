#ifndef THE_BASIC_C_HEADER_FILES
#define THE_BASIC_C_HEADER_FILES
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#endif

#ifndef GSL_HEADER_FILES
#define GSL_HEADER_FILES
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif

#include <GSL_stat.h>

#if defined CPGPLOT_REPRESENTATION
/* This header file contains the function prototypes
 * of my hand-made cpgplot-based plotting procedures
 */
#include <CPGPLOT_HEADER.h>
#endif


typedef struct Parameter_Modelinfo
{

#include <include.Parameter_Model.global.h>  
  
}Parameter_Model;

void ArgumentControl(int , char **); /* Only called from the main function 
						in main.c 
					     */

int numerical_Integration_Driver( Parameter_Model * ,
				  int , double * ,
				  double , double , 
                                  double * );
int function (double , const double [], double [], void *);
int jacobian (double , const double [], double *, double [], void * ); 
