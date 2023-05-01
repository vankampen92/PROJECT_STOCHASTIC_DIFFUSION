/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                            David Alonso, 2017 (c)                         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
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

#if defined NUMERICAL_RECIPES
#ifndef NUMERICAL_RECIPES_LIBRARY
#define NUMERICAL_RECIPES_LIBRARY
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Numerical Recipes Header Files:                                           */
#include <nr.h>
#include <nrutil.h>
#include <complex.h>
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#endif
#endif

#ifndef GSL_HEADER_FILES
#define GSL_HEADER_FILES
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#endif

#include <stat.h>

#if defined CPGPLOT_REPRESENTATION
/* This header file contains the function prototypes
 * of my hand-made cpgplot-based plotting procedures
 */
#include <CPGPLOT_HEADER.h>
#endif
