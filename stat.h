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

#define MAX(A,B) ((A)>(B)? A:B)
#define MIN(A,B) ((A)<(B)? A:B)

/* This file contains prototype definitions of mainly obsolete functions. Function code is on 
   stat.c. These utilities are useful, but they are mainly here to ensure back compatibility with 
   old code. They somehow are often used in conjuction with numerical recipes library, which 
   I have trying to avoid in my newer code in favor of the open source better alternatave 
   GNU Scientic Library (GSL). 
*/

/* Three common declation to make numerical recipes library to be built. They affect only a few
   functions from the library, which I have never used. 
*//* * * * * * * * * * * * * * * * * */
// long int idum;
// unsigned long ija[1];
// double sa[1];
/* * * * * * * * * * * * * * * * * */

typedef struct positioninfo
{
  double x;
  double y;
}position;

struct point
{
  float x;
  float y;
};

void Press_Key();
void Print_Press_Key(int Stat_Bool, int Key, char * Message);
void Print(int, char *);

/* When are two double variables approximately 
   equal within some level of accuracy ? */
int X_apx_Y( double , double , double );

void Apply_Negative_Control(double * y, int K); 

void Calculation_of_Local_Machine_MAX_min_allowed_float_numbers (float * xmin, float * xmax);

void machar(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	    int *iexp, int *minexp, int *maxexp, float *eps, float *epsneg,
	    float *xmin, float *xmax);

/* M E M O R Y    C O R R U P T I O N   T E S T * * */
void STAT_memory_Corrupcion_Check_Utility();

/* Some useful simple mathematical functions */
int potencia(int, int);
double lnGamma(int n);
double Euclidean_Distance_between_Two_Positions(position p_YX, position p_2);
float Euclidean_distance(float,float, float,float); 
double distance(struct point, struct point);
int CombiNumber(int k, int n);
double logCombiNumber(double x, int y);
double log_Combinatorial_Number(int N, int n);
double logProchrammer(double , int );
void Finding_AdjacenceMatrix_Active(int, int **, float **);
void Finding_AdjacenceMatrix_Passive(int, int **, float **);
void show_IntegerMatrix(int **, int);
void show_FloatMatrix(float **, int);
void show_a_view_FloatMatrix(float **, int, int, int, int);
void show_DoubleMatrix(double **, int, int);
void show_DoubleMatrix_Name(double **M, int Nx, int Ny, char * A);
void show_a_vector_in_a_row_Int(int *Vector, int No);
void show_a_vector_in_a_row_Float(float *Vector, int No);
void show_a_vector_in_a_row_Double(double *Vector, int No);
/* 
   Functions to sample data sets following particular statistical distributions 
*/
double log_normal_sampling(float, float);
float gasdev_01();
double SamplingBeta_Distribution(double, double);
int sample_NegBinom(double, double);
int sample_NegBinom_Truncate(double, double, int);
int sample_ZSM_Local(int, double, double);
double log_NegBinom_Prob(double q, double nu, int n);
double Log_NegBinom_Prob(double, double, int n);
double NegBinom_Prob(double, double, int);
double Log_Poisson(double, int);
double Binomial_Betai(int N, double p, int k);
double Binomial(int, double, int);
double log_Binom_Prob(int N, double p, int n);
double Combinatorial_Number(int N, int n);
double Poisson(double, int);
int random_int_drand48(int,int);
long int random_int_lrand48(int,int);
double Gamma_Density_Ditribution(float n, float g, double tau);

/* Functions in ./NumRecip/gamma_Double.c */
double gammQ_double(double a, double x);
double gammP_double(double a, double x);
double gcf_D(double a, double x);
double G_6_2_7(double a, double x);
double gser_D(double a, double x);
double G_6_2_5(double a, double x);

/* In ./NumRecip/gamma_Double.c  */
double Conditional_Gamma_Density_Distribution(double t, double n, double x);
/* 
   Computing simple statistical descriptors from data
*/
float Average(int [], float [], int); 
float Average_float_Vector(float *Pop, int Sp);
double Average_double_Vector(double * Pop, int Sp);
double Average_int_Vector(int *Pop, int Sp);
double average_squared_error_per_point(double * C, double * y, int N);
double average_absolute_error_per_point(double * C, double * y, int N);
double Variance_double_Vector(double *Pop, int Sp);
void Average_and_Variance_double_Vector(double *Pop, int Sp, double *Ave, double *Var);
double Correlation_double_Vector(double * y_0, double * y_1, int N);
void Average_Variance_and_Absolute_Deviation_double_Vector(double *Pop, int i_0, int Sp, double *Ave, double *Var, double *aDev);
float Variance_float_Vector(float *Pop, int Sp);
float Shannon(int [], int);
int Sum_iVector(int *, int );
double Sum_dVector(double *, int);
int Sum_Diversity(int *, int);
double Sum_Diversity_double(double *mVector, int Sp);
double mean_amplitude(float *X, int No);
double mean_squared_amplitude(float *X, int No);
void moving_Average_double_Vector(double * Pop, double * movPop, int Sp, int N);

/* 
   Hitograms and empirical distributions from real data 
*/
void abundance_distribution_int(int *, int, int *, int);
void abundance_distribution(int [], int, double [], int);
void abundance_distribution_2D(int [], double **, int);
void histogram_x_axes(double [],int, double,double);
void histogram(double [],int, double [],int, double,double);
void histograma(int, double [], int);
void revert(int, float []);
void rank_Abundance_0(int [], int, double []);
void rank_Abundance_Pop_YX_n(int [], int, double []);
void rank_Abundance_Ave(double [], int, double [], double []);
void rank_Abundace_Ave_Var(double [], double [], int , int);
void Ave_Var(double *, double *, int , int);
void accum_Diversity_Histogram(int *, int, double *, double *);
void accum_Histogram(double *, int, double *, double *);
void Norma(double [], int);
void Norma_Total_a_J(float *, int, int);
void probability_distribution_from_stochastic_realizations(double * y,
							   int Realizations,
							   double * p, int p_DIM );
void probability_distribution_from_stochastic_realizations_2D(double * x, double * y,
							      int Realizations,
							      double ** p,
							      int p_DIM_x, int p_DIM_y );
void Norma_2D_Nx_Ny(double **Relative_Abundance, int No_x, int No_y); 
void Norma_Total_0_J(double *, int);
void Norma_Total_0_J_Percentage(double *, int);
void Norma_2D(double **, int);
void Norm_P_n_m(int, double **P, double);
int Number_of_Octaves(int max_A, double * x, int MAX_NUMBER_of_OCTAVES_EVER);

/* 
   Simple operations with vectors and matrices
*/
void dvector_to_float(double *, float *, int);
void offsetting_float(float *, int);
void set_matrix_to_value_float (float **pN, int S, float value);
void set_matrix_to_value_float_ColumnsRows(float **, int, int, float);
void set_matrix_to_value_float_NxNy (float **, int, int, float);
void set_matrix_to_value_double_1_Nx_1_Ny (double **pN, int Nx, int Ny, double value);
void set_matrix_to_value_float_1_Nx_1_Ny (float **pN, int Nx, int Ny, float value);
void set_matrix_to_value_double(double **pN, int S, double value);
void set_matrix_to_value_double_ColumnsRows(double **, int, int, double);
void set_matrix_to_value_double_NxNy (double **, int, int, double);
void set_matrix_to_value_int_ColumnsRows(int **, int, int, int);
void set_to_value_double(double *pN, int S, double value);
void set_to_value_double_1_N (double *pN, int S, double value);
void set_to_value_float_1_N (float *pN, int S, float value);
void set_to_value_float(float *pN, int S, float value);
void set_to_value_int(int *pN, int S, int value);
int discrete_MAX_double(double *, int);
double vector_MAX_double(double *, int, int);
double vector_MAX_Index_double(double * , int , int , int * );
int vector_MAX_int(int *, int, int);
double vector_MIN_double(double *, int, int);
void herror(int, int, int, int);
void ordering_Vector_double(double *X, int n, int signe);
void ordering_Vector_2_double(double *X, double *Y, int n, int HOW);
int Finding_a_Value_int(const int * N, int S, int value,  int * i_th);

/* 
   Input --- Output
   Simple functions to save and read from disk 
   two vectors in a predefined, ordered file in a one, two column way, n
   column way, etc.
*/
void name_Ordered(char Prefix[], int File_Number, char Sufix[], char *File);
void label_Name( char * Name_OUT, char * Name_IN );
void fileName (char [], char [], char []);
void fitxer (char [], char [], int, char []);
void inttochar(int, char []);
void create_Title_MONTH_YEAR(char * Title, char * Pre_Title, int i, int Initial_Year, int Initial_Month);
int doubletochar(double, char *);
void Saving_to_File_int(char [], double [], int [], int, int);
void Saving_to_File_float_float(char [], float [], float [], int, int, int);
void Saving_to_File_float(char [], double [], float [], int, int);
void Saving_to_File_double(char [], double [], double [], int, int);
void Saving_to_File_double_1COL(char *Prefix, double Number[], int Points, int File_Number);
void Saving_to_File_double_1Col(char *, double *, int);
void Saving_to_File_double_3Col(char [], double [], double [], double [], int, int);
void Saving_to_File_double_4Col(char Prefix[], 
				double Time[], double y_1[], double y_2[], double y_3[], 
				int Points, int File_Number);
void Saving_to_File_double_Average_plusminus_VAR(char Prefix[], double Time[], double Number[], 
						 double Var[], int Points, int File_Number);
void Saving_to_File_distribution(char [], double [], int, int);
void Saving_to_File_Lattice_double(char *Prefix, double *x, double *y, double **Lattice, int Nx, int Ny, 
				   int File_Number);
void Saving_to_File_double_y(char * File_Name, double ** y_Data, int No_Rows, int No_Cols);
void save_a_Matrix_to_File_int(char[], int, int **, int, int);
void save_a_Matrix_to_File_float(char[], int, float **, int, int);
void save_a_TimeMatrix_to_File_float(char [],int, double [], float **,int,int);
void save_a_TimeMatrix_to_File_int(char [],int, double [], int **,int,int);
void Saving_floatVector(char *, float *, int);
void Saving_intVector(char*, int *, int);
void saving_Time_Two_Float_vectors(char *file, 
				   float *time, float *x, float *y, int nP);
void Reading_from_File_double_TwoCol(char *, double *, double *, int);
void Reading_from_File_double_nxy(char * File_Name, double * x_Data, double ** y_Data, int No, int n_XY);
void Reading_from_File_double_nxy_Counting_Rows(char * ,
						double * , double ** , int ,
						int , int * );
void Reading_from_File_double_y(char * , double ** , int , int );
void Reading_floatVector(char *, float *, int);
void Reading_DiversityData(char*, int *, int);
void Reading_DiversityData_double(char* File, int *Species, int Data);
void Reading_DiversityData_double_BIS(char* File, double *Species, int Data);
void Reading_DiversityData_ONE_COLUMN(char * , int *, int);
void Reading_Diversity_ONE_COLUMN(char *, int *, int Data);
void Reading_doubleVector(char* File, double *mVector, int Sp);
void Reading_floatVector(char* File, float *mVector, int Sp);
void Reading_intVector(char* File, int *mVector, int Sp);

/* Simple solution of two dimension systems of linear equations */
void linearEqSystem_2dim(double *p0, double *p1, double *p2, double *Sol);

