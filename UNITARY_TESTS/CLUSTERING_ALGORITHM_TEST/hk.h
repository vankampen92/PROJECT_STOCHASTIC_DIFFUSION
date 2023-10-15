/* Tobin Fricke's implementation of the Hoshen-Kopelman algorithm for
   cluster labeling (header file version)

   Copyright (c) September 9, 2000, by Tobin Fricke <tobin@splorg.org>
   Distributed under the terms of the GNU Public License.

   Modified 2002-03-09 Tobin Fricke
   Modified substantially 2004-04-21 by Tobin Fricke

   This program is written in the 1999 standard of the C language (C99).  Older C
   compilers will refuse to compile it.   You can use a C++ compiler, a C99 compiler,
   or you can modify this code to comply with a previous version of the C standard.
   The GCC compiler supports C99 as of version 3.0.  Compile this program with:

   gcc-3.0 -Wall -std=c99 -c hk.c

   http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
*/

/*  uf_find returns the canonical label for the equivalence class containing x */
int uf_find(int x);

/*  uf_union joins two equivalence classes and returns the canonical label of the resulting class. */
int uf_union(int x, int y);

/*  uf_make_set creates a new equivalence class and returns its label */
int uf_make_set(void);

/*  uf_intitialize sets up the data structures needed by the union-find implementation. */
void uf_initialize(int max_labels);

/*  uf_done frees the memory used by the union-find data structures */
void uf_done(void);

/* End Union-Find implementation */
#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

/* print_matrix prints out a matrix that is set up in the "pointer to pointers" scheme
   (aka, an array of arrays); this is incompatible with C's usual representation of 2D
   arrays, but allows for 2D arrays with dimensions determined at run-time */
void print_matrix(int **matrix, int m, int n);

/* Label the clusters in "matrix".  Return the total number of clusters found. */
int hoshen_kopelman(int **matrix, int m, int n);

/* This procedure checks to see that any occupied neighbors of an occupied site
   have the same label. */
void check_labelling(int **matrix, int m, int n);
