#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#endif

#if defined CPGPLOT_REPRESENTATION
/* This header file contains the function prototypes
 * of my hand-made cpgplot-based plotting procedures
 */
#include <CPGPLOT_HEADER.h>
#endif

#include <BP_StrainsProfiles.h>
