#include "stat.h"

void Apply_Negative_Control(double * y, int K)
{
  /* Some times negative values make no sense...
     Let us get rid of them */
  int i;
  for(i=0; i<K; i++) y[i] = MAX(y[i], 0.0);
}


int X_apx_Y( double X, double Y, double ACU)
{
  /* When are two double variables approximately
     equal within some level of accuracy ? */

  int R;

  if ( fabs(X - Y) < ACU ) R = 1; /* The two number are appox the same within ACU */
  else                     R = 0; /* The two number differ in more than ACU       */

  return(R);
}

void Press_Key()
{
  printf("Press any key to continue...\n");
  getchar();
}

void Print_Press_Key(int Stat_Bool, int Key, char * Extra_Message)
{
  if(Stat_Bool == 1) printf("%s\n", Extra_Message);
  if(Key       == 1) {
    printf("Press any key to continue... or Ctr-C to kill it\n");
    getchar();
  }
}

void Print(int Stat_Bool, char *Message)
{
  if(Stat_Bool == 1) printf("%s\n",Message);
}

void herror(int flag, int n0, int Pop0, int iter)
{
  if(n0 != Pop0) {
      printf("%d Error... n0 != Pop[0], your are getting into serious trouble!!!\n",
	     iter);
      printf("n0 = %d, but Pop[0] = %d\tFlag = %d\n",n0,Pop0,flag);
  }
}

void STAT_memory_Corrupcion_Check_Utility()
{
  /* M E M O R Y    C O R R U P T I O N   T E S T * * */

  /* Usually, when memory corruption has occurred without
     a crash, and so without notice... further allocation-
     de-allocation requests make the running program
     totally crash (segmentation fault) because memmory
     handling has been corrupted. By calling this
     funciton all over your code, you will be able to track
     down when the unadverted corruption actually occurred.

     You may call this function by writing an only line of code:

     STAT_memory_Corrupcion_Check_Utility();
  */
  /* This function is part of the my stat.c library.
  */
  int * merda; static int quanta_merda = 0;
  printf("/* * * * * * * * * * * * * * * * * * * * * * * * */\n");
  printf("/* M E M O R Y   C O R R U P T I O N   T E S T   */\n");
  printf("merda allocation trial = %d\n", quanta_merda);
  merda = (int *)malloc( 100 * sizeof(int) );
  printf("merda allocation succeeded = %d\n", quanta_merda);
  printf("merda DE-allocation trial = %d\n", quanta_merda);
  free (merda) ;
  printf("merda de-allocation succeeded = %d\n", quanta_merda);
  printf("/* M E M O R Y   I S   I N   G O O D   S H A P E */\n");
  printf("/* * * * * * * * * * * * * * * * * * * * * * * * */\n");

  quanta_merda++;
  /* * * * * * * * * * * * * * * * * * * * * * * * * */
}

void Calculation_of_Local_Machine_MAX_min_allowed_float_numbers (float * xmin, float * xmax)
{
  /* Evaluation of the largest (xmax) and the tiniest (xmin)
     numbers my machine can handle
  */
    int ibeta; int it; int irnd; int ngrd; int machep; int negep;
    int iexp; int minexp; int maxexp; float eps; float epsneg;

    machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep,
           &iexp, &minexp, &maxexp, &eps, &epsneg,
           xmin, xmax);
}

/* B E G I N : Some functions to calculate easy statistics */
void Finding_AdjacenceMatrix_Passive(int Sp, int **A, float **M)
{
  /* Finding the adjacent matrix... */
  int i,j,ii;

  for (i=0; i<=Sp; i++){
    ii = 0;
    for (j=0; j<=Sp; j++)
      if(M[i][j] >= M[j][i])
	if( MAX(M[i][j], M[j][i]) > 0){        /* j-th species feeds on i-th species */
	  A[i][ii] = j;
	  ii++;
	}
    for(j = ii; j<Sp; j++) A[i][j] = 0;
    A[i][Sp] = ii ;                 /* Number of links of species i */
  }
}

void Finding_AdjacenceMatrix_Active(int Sp, int **A, float **M)
{
  /* Finding the adjacent matrix... */
  int i,j,ii;

  for (i=0; i<=Sp; i++){
    ii = 0;
    for (j=0; j<=Sp; j++)
      if(M[i][j] <= M[j][i])
	if( MAX(M[i][j], M[j][i]) > 0){        /* j-th species is food for i-th species */
	  A[i][ii] = j;
	  ii++;
	}
    for(j = ii; j<Sp; j++) A[i][j] = 0;
    A[i][Sp] = ii ;                 /* Number of links of species i */
  }
}

void show_IntegerMatrix(int **M, int Sp)
{
  int i,j, Nx,Ny, n_sys;

  /* n_sys = system("clear"); */

  //Nx = Sp+1 > 27 ? 27:Sp+1;
  //Ny = Sp+1 > 17 ? 17:Sp+1;

  Nx = Sp > 27 ? 27:Sp;
  Ny = Sp > 17 ? 17:Sp;
  printf("Printing adjacent matrix...\n");

  printf("* :\t");
  for(i=0; i<Nx; i++) printf("%3d  ",i);
  printf("\n\n");

  for(i=0; i<Nx; i++){
    printf("%d :\t",i);
    for(j=0; j<Ny; j++) printf("%3d  ",M[i][j]);
    printf("\n");
  }
  printf("\n\n");
  /* getchar(); */
}

void show_FloatMatrix(float **M, int Sp)
{
  int i,j, Nx,Ny, n_sys;

  /*n_sys = system("clear");*/
  printf("Printing colonization matrix...\n");

  Nx = Sp+1 > 27 ? 27:Sp+1;
  Ny = Sp+1 > 17 ? 17:Sp+1;

  printf("* :\t");
  for(i=0; i<Nx; i++) printf("%5d  ",i);
  printf("\n\n");

  for(i=0; i<Nx; i++){
    printf("%d :\t",i);
    for(j=0; j<Ny; j++) printf("%4.3f  ",M[i][j]);
    printf("\n");
  }
  printf("\n\n");
  /* getchar(); */
}

void show_a_view_FloatMatrix(float **M, int Nx_0, int Nx, int Ny_0, int Ny)
{
  int i,j;

  printf("* :\t");
  for(i=Nx_0; i<=Nx; i++) printf("%5d  ",i);
  printf("\n\n");

  for(i=Nx_0; i<=Nx; i++){
    printf("%d :\t",i);
    for(j=Ny_0; j<=Ny; j++) printf("%4.3f  ",M[i][j]);
    printf("\n");
  }
  printf("\n\n");
  /* getchar(); */
}

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

void show_DoubleMatrix_Name(double **M, int Nx, int Ny, char * A)
{
  int i,j, n_sys;

  /*n_sys = system("clear");*/
  //printf("Printing matrix...\n");

  printf("%s :\t", A);
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

void show_a_vector_in_a_row_Int(int *Vector, int No)
{
  int i;
  for(i=0; i<No; i++)
   if(i<No-1)
      printf("%d, ", Vector[i]);
    else
      printf("%d.\n", Vector[i]);
}

void show_a_vector_in_a_row_Float(float *Vector, int No)
{
  int i;
  for(i=0; i<No; i++)
   if(i<No-1)
      printf("%f, ", Vector[i]);
    else
      printf("%f.\n", Vector[i]);
}

void show_a_vector_in_a_row_Double(double *Vector, int No)
{
  int i;
  for(i=0; i<No; i++)
   if(i<No-1)
      printf("%g, ", Vector[i]);
    else
      printf("%g.\n", Vector[i]);
}

int potencia(int a, int n)
{
  /* An integer a raised to the power of an integer n */
  int i;
  int R;

  R = 1;
  for(i=0; i<n; i++)
    R *= a;

  return R;
}

int CombiNumber(int n, int k)
{
  int a;
  double logof_a, x; 

  if(n < k){
    printf("No Combinatorial Number Exists...\n");
    printf("n = %d\tk = %d\n",n, k);
    printf("n is less than k!!!!\n");
    exit(0);
  }
  else if(k == n)
    a = 1;
  else{
    x       = (double)n;
    logof_a = logCombiNumber(x, k);     
  }

  a = (int)exp(logof_a);

  return(a);
}

double logCombiNumber(double x, int y)
{
  double a;
  int k;

  if(x < y){
    printf("No Combinatorial Number Exists...\n");
    printf("x = %g\ty = %d\n",x,y);
    printf("x is less than y!!!!\n");
    exit(0);
  }
  else if((int)x == y)
    a = 0.;
  else{
    a = 0.;
    for(k=0; k<y; k++)
      a += log((x-(double)k)/(y-(double)k));
  }

  return(a);
}

double logProchrammer(double I, int n)
{
  double a;
  int k;

  if(I <= 0 || n <0){
    printf("No Prochrammer number exists...\n");
    printf("I = %g\tn = %d\n",I,n);
    printf("I is less than 0!!!!\n");
    exit(0);
  }
  else if(n == 0)
    a = 0.;
  else{
    a = 0.;
    for(k=0; k<n; k++)
      a += log((I+(double)k));
  }

  return(a);
}

#if defined NUMERICAL_RECIPES
double lnGamma(int n)
/*
   The factorial of the integer $n-1$ approximated by the gamma function when it is greater
   than 10 and computed exactly when lower. The output is given as the logarithm of that
   factorial
*/
{
  double f;
  int m;
  /* assert(n>0); */

  if(n < 10){
    f = 1.;
    m = n-1;
    while(m>1){
      /* printf("%d\t%d\t%g\n",n,m,f); */
      f *= (double)m--;
    }
    f = log(f);
  }
  else
    f = (double)gammln((float)n);

  return f;
}
#endif

double distance(struct point A, struct point B)
{
  double d;
  d = sqrt((A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y));
  return d;
}

double Euclidean_Distance_between_Two_Positions(position p_YX, position p_2)
{
  double Distance;
  Distance = sqrt((p_YX.x - p_2.x)*(p_YX.x - p_2.x) + (p_YX.y - p_2.y)*(p_YX.y - p_2.y) );
  return Distance;
}

float Euclidean_distance(float x1,float y1, float x2,float y2)
{
  float Distance;
  Distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  return Distance;
}

double mean_amplitude(float *X, int No)
{
  /* This is the mean value of the series X[] */
  int i;
  double mean;

  mean = 0.;
  for(i=1; i<=No; i++) mean += X[i];

  mean /= (double)No;

  return mean;
}

double mean_squared_amplitude(float *X, int No)
{
  /* This is the mean of squared values of the
     series X[]. This is not the variance */
  int i;
  double mean;

  mean = 0.;
  for(i=1; i<=No; i++) mean += X[i]*X[i];

  mean /= (double)No;

  return mean;
}

float Average_float_Vector(float *Pop, int Sp)
{
    /* Notice that the vector starts at i=1 */
    int i;
    float Sum, Ave;

    Sum = 0.;
    for (i=1; i<= Sp; i++) Sum += Pop[i];

    Ave = Sum/(float)Sp;

    return Ave;
}

double Average_int_Vector(int *Pop, int Sp)
{
    /* Notice that the vector starts at i=1 */
    int i;
    double Sum;
    double Ave;

    Sum = 0.;
    for (i=1; i<= Sp; i++) Sum += (double)Pop[i];

    Ave = Sum/(double)Sp;

    return Ave;
}

void moving_Average_double_Vector(double * Pop, double * movPop, int Sp, int N)
{
    /* Notice that the vector starts here at i=0                */
    /* N, order of the moving average                           */
    /* Achtung: the order must be an even (2, 4, ...) number!!! */
    int i,j;
    double Sum, Ave;

    if(N%2 == 0){
      for (i=N/2; i< Sp-N/2; i++){

	Sum = 0.;
	for (j = i-N/2; j < i+N/2; j++) Sum += Pop[j];

	Ave = Sum/(double)N;

	movPop[i] = Ave;
      }
    }
    else{
      printf(" In this implementation (moving_Average_double_Vector() in stat.c)\n");
      printf(" The order of the moving average most be an\n");
      printf(" even number. This program will self-desctruct now!!\n");
      exit(0);
    }

    /* Completing the ends of the vector with movPop[N/2] and movPop[Sp-N/2-1] */

    for(i=0; i<N/2; i++) movPop[i] = movPop[N/2];

    for(i=Sp-N/2; i< Sp; i++) movPop[i] = movPop[Sp - N/2 - 1];
}

double Average_double_Vector(double *Pop, int Sp)
{
    /* Notice that the vector starts here at i=0*/
    int i;
    double Sum, Ave;

    Sum = 0.;
    for (i=0; i < Sp; i++) Sum += Pop[i];

    Ave = Sum/(double)Sp;

    return Ave;
}

double average_squared_error_per_point(double * C, double * y, int N)
{
  int i;
  double x;

  x = 0.0;

  for (i=0; i < N; i++){

    x += (C[i] - y[i])*(C[i] - y[i]);

  }

  x /= (double)N;

  return(x);
}

double average_absolute_error_per_point(double * C, double * y, int N)
{
  int i;
  double x;

  x = 0.0;

  for (i=0; i < N; i++){

    x += fabs(C[i] - y[i]);

  }

  x /= (double)N;

  return(x);
}

void Average_and_Variance_double_Vector(double *Pop, int Sp, double *Ave, double *Var)
{
  /* Notice that the vector starts at i=0 (Como Dios ordena y manda) */
    int i;
    double Sum;

    (*Ave) = 0.; (*Var) = 0.;

    for (i = 0; i < Sp; i++){
      (*Ave) += Pop[i];
      (*Var) += Pop[i]*Pop[i];
    }

    if( Sp > 0 )  {
      (*Ave) = (*Ave)/(double)Sp;
      (*Var) = (*Var)/(double)Sp;
    }

    (*Var) = (*Var) - (*Ave)*(*Ave);
}

double Correlation_double_Vector(double * y_0, double * y_1, int N)
{
  int i;
  double corr;

  double AVE_0, VAR_0;
  double AVE_1, VAR_1;

  Average_and_Variance_double_Vector(y_0, N, &AVE_0, &VAR_0 );
  Average_and_Variance_double_Vector(y_1, N, &AVE_1, &VAR_1 );

  corr = 0.0;
  for(i = 0; i<N; i++) {
    corr += (y_0[i] - AVE_0) * (y_1[i] - AVE_1) / (double)N;
  }

  corr = corr / sqrt(VAR_0 * VAR_1);

  return (corr);
}

void Average_Variance_and_Absolute_Deviation_double_Vector(double *Pop, int i_0, int Sp, double *Ave, double *Var, double *aDev)
{
  /* Notice that the vector starts at i=0 (Como Dios ordena y manda) */
  int i;
  double Sum;

    (*Ave) = 0.0; (*Var) = 0.0; (*aDev) = 0.0;

    for (i = i_0; i < Sp; i++){
      (*Ave) += Pop[i];
      (*Var) += Pop[i]*Pop[i];
    }
    (*Ave) = (*Ave)/(double)Sp;

    (*Var) = (*Var)/(double)Sp;
    (*Var) = (*Var) - (*Ave)*(*Ave);

    if( (*Var) < 0) {
      for(i=i_0; i<Sp; i++) printf("%g, ", Pop[i]);
    }

    assert( Sp > 0);
    assert( (*Var) >= 0 );

    for (i = i_0; i < Sp; i++){
      (*aDev) += fabs( Pop[i] - (*Ave) );
    }
    (*aDev) = (*aDev)/(double)Sp;
}

double Variance_double_Vector(double *Pop, int Sp)
{
    /* Notice that the vector starts at i=0 */
    int i;
    double Sum, Ave, Var;

    Ave = 0.; Var = 0.;
    for (i = 0; i < Sp; i++){
      Ave += Pop[i];
      Var += Pop[i]*Pop[i];
    }
    Ave = Ave/(double)Sp;
    Var = Var/(double)Sp;

    Var = Var - Ave*Ave;

    return (Var);
}

float Variance_float_Vector(float *Pop, int Sp)
{
    /* Notice that the vector starts at i=1 */
    int i;
    float Sum, Ave, Var;

    Ave = 0.; Var = 0.;
    for (i=1; i<= Sp; i++){
      Ave += Pop[i];
      Var += Pop[i]*Pop[i];
    }
    Ave = Ave/(float)Sp;
    Var = Var/(float)Sp;
    Var = sqrt(Var - Ave*Ave);

    return Var;
}

float Average(int Pop[], float C[], int Sp)
{
  int i, kount;
  float x;

  kount = 0;
  for (i=1; i<= Sp; i++) kount += Pop[i];
  x = 0.;
  for (i=1; i<= Sp; i++)
    x += Pop[i]/(float)kount * C[i];
  return x;
}

/* Setting Values */
void set_to_value_double (double *pN, int S, double value)
{
  int i;
  for (i=0; i<S; i++) *pN++ = value;
}

void set_to_value_double_1_N (double *pN, int S, double value)
{
  int i;
  for (i=1; i<=S; i++) *pN++ = value;
}

void set_to_value_float_1_N (float *pN, int S, float value)
{
  int i;
  for (i=1; i<=S; i++) *pN++ = value;
}

void set_matrix_to_value_double (double **pN, int S, double value)
{
  int i,j;
  for (j=0; j<=S; j++)
    for (i=0; i<=S; i++) pN[i][j]= value;
}

void set_matrix_to_value_double_NxNy (double **pN, int Nx, int Ny, double value)
{
  int i,j;
  for (j=0; j<Ny; j++)
    for (i=0; i<Nx; i++) pN[i][j]= value;
}

void set_matrix_to_value_double_1_Nx_1_Ny (double **pN, int Nx, int Ny, double value)
{
  int i,j;
  for (j=1; j<=Ny; j++)
    for (i=1; i<=Nx; i++) pN[i][j]= value;
}

void set_matrix_to_value_float_1_Nx_1_Ny (float **pN, int Nx, int Ny, float value)
{
  int i,j;
  for (j=1; j<=Ny; j++)
    for (i=1; i<=Nx; i++) pN[i][j]= value;
}

void set_matrix_to_value_double_ColumnsRows(double **pN, int Rows, int Columns, double value)
{
  int i,j;
  for (j=0; j<Columns; j++)
    for (i=0; i<Rows; i++) pN[i][j]= value;
}

void set_matrix_to_value_float(float **pN, int S, float value)
{
  int i,j;
  for (j=0; j<=S; j++)
    for (i=0; i<=S; i++) pN[i][j]= value;
}

void set_matrix_to_value_float_NxNy (float **pN, int Nx, int Ny, float value)
{
  int i,j;
  for (j=0; j<Ny; j++)
    for (i=0; i<Nx; i++) pN[i][j]= value;
}

void set_matrix_to_value_float_ColumnsRows(float **pN, int Rows, int Columns, float value)
{
  int i,j;
  for (j=0; j < Columns; j++)
    for (i=0; i < Rows; i++) pN[i][j]= value;
}

void set_matrix_to_value_int_ColumnsRows(int **pN, int Rows, int Columns, int value)
{
  int i,j;
  for (j=0; j < Columns; j++)
    for (i=0; i < Rows; i++) pN[i][j]= value;
}

void set_to_value_float (float *pN, int S, float value)
{
  int i;
  for (i=0; i<S; i++) *pN++ = value;
}

void set_to_value_int (int *pN, int S, int value)
{
  int i;
  for (i=0; i<S; i++) *pN++ = value;
}

void dvector_to_float(double *y, float *f_y, int n)
{
  int i;

  double * dummy_y;
  float  * dummy_f_y;

  dummy_y   = y;
  dummy_f_y = f_y;

  for(i=0; i<n; i++) {
    *dummy_f_y++ = *dummy_y++;
    /* f_y[i] = (float)y[i]; */
    /* if(f_y[i] == 0.0) {
         printf("vector(%d) = %f\t%g\n",i,f_y,y);
         exit(0);
       }
    */
  }
}

void offsetting_float(float *y, int n)
{
  int i;
  for(i=n; i>0; i--) y[i] = y[i-1];
}

/*
   Functions to sample data sets following particular statistical distributions
*/
double SamplingBeta_Distribution(double u, double beta)
{
  /* This is the Beta distribution of parameters (u, beta) just when u is equal to 1. */
  double x, L;

  assert(u == 1.);
  L = 1./beta * log(1. - drand48());
  x = 1. - exp(L);

  return x;
}

long int random_int_lrand48(int a1, int a2)
{
  /* A random integer is generated from n = a1 to n = a2, both extrems permitted */
  int n;
  long int n_random;

  n = a2 - a1 + 1;
  n_random = a1 + lrand48()%n;

  return (n_random);
}

int random_int_drand48(int a1, int a2)
{
  /* A random integer is generated from n = a1 to n = a2, both extrems permitted */
  int n;
  int n_random;
  int stat_Bool;

  n = a2 - a1 + 1;
  stat_Bool = 0;
  while (stat_Bool == 0){

    n_random = (int)ceil( (double)a1 + drand48() * (double)n );

    if(n_random <= a2) stat_Bool = 1;
  }

  return (n_random);
}

int sample_NegBinom(double lambda, double nu)
{
  /* lambda is the MEAN and nu is the CLUMPLING parameter  */
  int i;
  float x;
  int deviate;
  long idum = -6509;
  double q,F, Sum;

  if(nu > 0.){
    q = lambda/nu;
    x = drand48();
    i = 0;
    Sum = 0;
    while(Sum < x){
      Sum += NegBinom_Prob(q, nu, i++);
    }
    deviate = i-1;

    return (deviate);
  }
  else
    return (0);
}

int sample_NegBinom_Truncate(double lambda, double nu, int N_max)
{
  int i;
  float x;
  int deviate;
  long idum = -6509;
  double q,F, Sum;

  q = lambda/nu;
  x = drand48();
  i = 0;
  Sum = 0;
  while((Sum < x) && (i < N_max)){
    Sum += NegBinom_Prob(q, nu, i++);
  }
  deviate = i-1;

  return (deviate);
}

double Gamma_Density_Ditribution(float n, float g, double tau)
{
  double logGamma, p;

  logGamma = n*log((double)g*(double)n)-gammln(n)+(n-1)*log(tau)-g*n*tau;

  p = exp(logGamma);
  return (p);
}

double Poisson(double lambda, int n)
{
  double Pn;

  Pn = -lambda + (double)n * log(lambda) - gammln((float)n+1.);
  Pn = exp(Pn);
  return (Pn);
}

double Log_Poisson(double lambda, int n)
{
  double Pn;

  Pn = -lambda + (double)n * log(lambda) - gammln((float)n+1.);

  return (Pn);
}

#if defined NUMERICAL_RECIPES
double Binomial_Betai(int N, double p, int k)
{
  /* This is calculating by substrancting two accumulated binomials */
  double Bin;
  double P_k, P_k_1;

  P_k = betai(k, N-k+1,(float)p); P_k_1 = betai(k, N- k + 2,(float)p);

  Bin = P_k - P_k_1;

  return(Bin);
}
#endif

double Combinatorial_Number(int N, int n)
{
  int i;
  double logCombinatorial_Number;
  double number;

  logCombinatorial_Number = 0.;
  for(i=0; i<n; i++){
    logCombinatorial_Number += log((double)(N-i)) - log((double)(n-i));
  }
  number = exp(logCombinatorial_Number);

  return(number);
}

double log_Combinatorial_Number(int N, int n)
{
  int i;
  double logCombinatorial_Number;

  logCombinatorial_Number = 0.;
  for(i=0; i<n; i++){
    logCombinatorial_Number += log((double)(N-i)) - log((double)(n-i));
  }

  return(logCombinatorial_Number);
}

double Binomial(int N, double p, int n)
{
  /* This is calculating by using combinatorial numbers */

  double Bin;
  double logCombinatorial_Number;
  int i;

  logCombinatorial_Number = 0.;
  for(i=0; i<n; i++){
    logCombinatorial_Number += log((double)(N-i)) - log((double)(n-i));
  }

  Bin = logCombinatorial_Number + (double)n*log(p) + (double)(N-n)*log(1.-p);

  Bin = exp(Bin);

  return(Bin);
}

double log_Binom_Prob(int N, double p, int n)
{
  /* This is calculating by using combinatorial numbers */

  double Bin;
  double logCombinatorial_Number;
  int i;

  logCombinatorial_Number = 0.;
  for(i=0; i<n; i++){
    logCombinatorial_Number += log((double)(N-i)) - log((double)(n-i));
  }

  Bin = logCombinatorial_Number + (double)n*log(p) + (double)(N-n)*log(1.-p);

  return(Bin);
}

double NegBinom_Prob(double q, double nu, int n)
{
  int i;
  double value;
  double logFac;

  logFac = 0.;
  for(i=n; i>1; i--) logFac += log((double)i);

  value = gammln((float)nu+(float)n) - gammln((float)nu) - logFac + n*log(q) - (nu+n)*log(1+q);

  value = exp(value);
  return(value);
}

double Log_NegBinom_Prob(double q, double nu, int n)
{
  double P;
  double x;

  x = (double)n + nu - 1.;
  P = logCombiNumber(x,n) + (double)n*log(q) - ((double)n+nu)*log(1.+ q);

  return(P);
}

double log_NegBinom_Prob(double q, double nu, int n)
{
  int i;
  double value;
  double logFac;

  logFac = 0.;
  for(i=n; i>1; i--) logFac += log((double)i);

  value = gammln((float)nu+(float)n) - gammln((float)nu) - logFac + n*log(q) - (nu+n)*log(1+q);

  return(value);
}

double log_normal_sampling(float Mean, float Dvst)
{
  double x;
  float y;
  /* long idum = 13456; ---> idum must be a Global Variable */
  long idum = 13456;

  y = gasdev_01();
  x = (double)Mean * exp(Dvst*y);
  return x;
}

float gasdev_01()
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	//	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*drand48()-1.0;
			v2=2.0*drand48()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

float Shannon(int Abundance[], int Sp)
{
  int i,sum;
  float H,p;
  int *pA;

  pA = &Abundance[1];
  sum = 0;
  for(i=1;i<=Sp;i++) sum += *pA++;
  H = 0.;
  for(i=1;i<=Sp;i++)
    if(Abundance[i]>0){
      p = (float)Abundance[i]/(float)sum;
      H += -p*log(p);
    }
  return H;
}

/*
   Hitograms and empirical distributions from real data
*/
void revert(int n, float V[])
{
  int i;
  float x;
  for (i=0; i < n/2; i++){
    x = V[i];
    V[i] = V[n-i-1];
    V[n-i-1] = x;
  }
}

void abundance_distribution_non_optimized(int Pop[], int Sp, double Relative_Abundance[], int No)
{
  int i,j;
  for(i=1; i<=Sp; i++){
    if(Pop[i] > No)
      Relative_Abundance[No+1]++;
    else
      for(j=0; j<=No; j++){
	if(Pop[i] == j)
	  Relative_Abundance[j]++;
      }
  }
}

void abundance_distribution_int(int Pop[], int Sp, int * Relative_Abundance, int No)
{
  int i;
  for(i=0; i < Sp; i++){
    if(Pop[i] > No)
      Relative_Abundance[No+1]++;
    else /* if(Pop[i]>0) */
      Relative_Abundance[Pop[i]]++;
  }
}

void probability_distribution_from_stochastic_realizations(double * y, int Realizations,
							   double * p, int p_DIM )
{
  /* Input:
     . y is an array with a number of 'Realizations' of a stochastic variable that takes
     integer values.
     . p_DIM is the max dimension of the p array

     Output:
     . p is the distribution (normalized to one) with the relative appearance of those
     integer values in the array y.

     For this algorithm to work, p_DIM > MAX_i ( y[i] ), and p should have been initialized
     with calloc() to make sure all its entries are zero.
  */

  int i, Pop;
  for(i=0; i < Realizations; i++){
    Pop = (int)y[i];
    if(Pop < p_DIM) p[Pop]++;
    else {
      printf(" The dimension of the distribution vector, p, is too low (p_DIM = %d)\n", p_DIM);
      printf(" It is lower than required. p_DIM should be larger or equal to MAX_i ( y[i] )\n");
      printf(" Population at Realization %d-th: Pop = %d\n", i, Pop);
      printf(" This p vector should be allocated with a higher dimension.\n");
      printf(" The program will exit\n");
      exit(0);
    }
  }

  Norma_Total_0_J(p, p_DIM-1);
}

void probability_distribution_from_stochastic_realizations_2D(double * x, double * y,
							      int Realizations,
							      double ** p,
							      int p_DIM_x, int p_DIM_y )
{
  /* Input:
     . x is an array with a number of 'Realizations' of the first stochastic variable that takes
     integer values.
     . x is an array with a number of 'Realizations' of the first stochastic variable that takes
     . p_DIM_x is the max dimension of the p array along the first dimension
     . p_DIM_y is the max dimension of the p array along the second dimension

     To be clear, the result of the i-th replicate was

     (x[i], y[i])

     and so on and so forth from i=0 to i=Realizations-1

     Output:
     . p is the distribution (normalized to one) with the relative appearance of those
     integer values in the data (x, y).

     For this algorithm to work, p_DIM_x > MAX_i ( x[i] ), and p_DIM_y > MAX_i ( y[i] ),
     but 'p' should have been initialized with calloc() to make sure all its entries are
     zero.
  */

  int i, Pop_x, Pop_y;
  for(i=0; i < Realizations; i++){
    Pop_x = (int)x[i];
    Pop_y = (int)y[i];

    if(Pop_x < p_DIM_x && Pop_y < p_DIM_y )  p[Pop_x][Pop_y]++;
    else {
      printf(" The dimension of the distribution matrix, p, is too low (p_DIM_x = %d)\n", p_DIM_x);
      printf(" The dimension of the distribution matrix, p, is too low (p_DIM_y = %d)\n", p_DIM_y);
      printf(" It is lower than required\n");
      printf(" p_DIM_x should be larger or equal to MAX_i ( x[i] ) and \n");
      printf(" p_DIM_y should be larger or equal to MAX_i ( y[i] ).\n");
      printf(" This p matrix should be allocated with higher dimensions.\n");
      printf(" The program will exit\n");
      exit(0);
    }
  }

  Norma_2D_Nx_Ny(p, p_DIM_x, p_DIM_y );
}

void Norma_2D_Nx_Ny(double **Relative_Abundance, int No_x, int No_y)
{
  int i,j;
  double summ;

  summ = 0.;
  for(i=0; i< No_x; i++)
    for(j=0; j< No_y; j++)
      summ += Relative_Abundance[i][j];

  for(i=0; i< No_x; i++)
    for(j=0; j< No_y; j++)
      Relative_Abundance[i][j] /= summ;
}



void abundance_distribution(int Pop[], int Sp, double Relative_Abundance[], int No)
{
  int i;
  for(i=1; i<=Sp; i++){
    if(Pop[i] > No)
      Relative_Abundance[No+1]++;
    else /* if(Pop[i]>0) */
      Relative_Abundance[Pop[i]]++;
  }
}

void abundance_distribution_2D(int Pop[], double **Relative_Abundance, int No)
{
  int i,j;
  int m, n;

  m = Pop[0] + Pop[1];
  n = Pop[1];

  for(i=0; i<=No; i++)
    for(j=0; j<=No; j++)
      if(n == i && m == j)
	Relative_Abundance[n][m]++;

  if(n > No || m > No)
    Relative_Abundance[No+1][No+1]++;
}

void histograma(int i, double Relative_Probability[], int No)
{
  if(i <= No)
    Relative_Probability[i]++;
  else if(i > No)
    Relative_Probability[No+1]++;
}

void histogram_x_axes(double X[], int noBar, double x_i, double x_s)
{
  int i;
  double step;

  step = (x_s - x_i)/(double)noBar;
  for (i=0; i<noBar; i++)
    X[i] = x_i + (step*i + step*(i+1))/2.;
  X[noBar] = x_i + step*noBar;
}

void histogram(double Time[], int no, double Pn[], int B, double x_i, double x_s)
{
  int i,j, logic;
  double step;

  step = (x_s - x_i)/(double)B;
  for(i=0; i<no; i++){
    logic = 0;
    j = 0;
    while(logic == 0 && j < B){
      if(Time[i] > (x_i + step*j) && Time[i] <= (x_i + step*(j+1)) ){
	      logic = 1;
	      Pn[j]++;
      }
      j++;
    }
    if(j == B) Pn[B]++;
  }
}

void Norma_Total_a_J(float *Relative_Abundance, int a, int J)
{
  /* Normalizing a vector defined form i=a to i=J
     Usually $a$ is either a = 0 or a = 1
  */
  int i;
  double summ;
  float *Dumm;
  summ = 0.; Dumm = &Relative_Abundance[a];
  for(i=a; i<= J; i++) summ += *Dumm++;

  if(summ > 0.){
    Dumm = &Relative_Abundance[a];
    for(i=a; i<= J; i++) *Dumm++ /= summ;
  }
  else{
    printf("ERROR... Summ = 0 --> Non-normalizable probability\n");
    exit(0);
  }
}

void Norma_Total_0_J(double *Relative_Abundance, int J)
{
  int i;
  double summ;
  double *Dumm;
  summ = 0.; Dumm = Relative_Abundance;

  for(i=0; i<= J; i++) summ += *Dumm++;
  if(summ > 0.0){
    Dumm = Relative_Abundance;
    for(i=0; i<= J; i++) *Dumm++ /= summ;
  }
  else{
    printf("ERROR... Summ = 0 --> Non-normalizable probability\n");
    exit(0);
  }
}

void Norma_Total_0_J_Percentage(double *Relative_Abundance, int J)
{
  int i;
  double summ;
  double *Dumm;

  summ = 0.; Dumm = Relative_Abundance;
  for(i=0; i<= J; i++) summ += *Dumm++;
  if(summ > 0.){
    Dumm = Relative_Abundance;
    for(i=0; i<= J; i++){
      *Dumm *= 100.;
      *Dumm++ /= summ;
    }
  }
  else{
    printf("ERROR... Summ = 0 --> Non-normalizable probability\n");
    exit(0);
  }

}

void Norma(double Relative_Abundance[], int No)
{
  int i;
  double summ;
  summ = 0.;
  for(i=0; i<= No+1; i++) summ += Relative_Abundance[i];
  if(summ > 0.0)
    for(i=0; i<= No+1; i++) Relative_Abundance[i] /= summ;
  else{
    printf("ERROR... Summ = 0 --> Non-normalizable probability\n");
    exit(0);
  }
}

void Norm_P_n_m(int ND, double **P, double Norm)
{
  int j;
  double *pP;

  ND = ND*ND;
  /* Normalising probability... */
  pP= P[0];
  for(j=0; j<=ND; j++,pP++)
      (*pP) /= Norm;
}

void Norma_2D(double **Relative_Abundance, int No)
{
  int i,j;
  double summ;
  summ = 0.;
  for(j=0; j<= No; j++)
    for(i=0; i<= No; i++) summ += Relative_Abundance[i][j];
  summ += Relative_Abundance[No+1][No+1];

  for(j=0; j<= No; j++)
    for(i=0; i<= No; i++) Relative_Abundance[i][j] /= summ;
}

#if defined NUMERICAL_RECIPES
void rank_Abundance_Pop_YX_n(int Pop[], int n, double Rank[])
{
  /* One thousand is the maximum number of species in the system */
  int i,m;
  float *Pop_O;

  Pop_O = vector(1,n);
  m = 0;
  for (i=1; i<=n; i++)
    if(Pop[i] > 0)
      Pop_O[++m] = -Pop[i];
  sort(m,Pop_O);  /* revert(n+1,Pop_O); */
  for (i=0; i<m; i++) Rank[i] += -(double)Pop_O[i+1];
  //for (i=m; i<n; i++) Rank[i] += 0.;
  /* Rank[0] is storing the maximum value of the n-dimensional vector
     (Pop[1],...,Pop[n]), where Rank[i=n] would be actually always equal to zero,
     if Rank[i] were defined from i=0 to i=n.
  */
  free_vector(Pop_O, 1,n);
}

void rank_Abundance_0(int Pop[], int n, double Rank[])
{
  /* Notice that Pop[] must be non-negative vector */
  int i,m;
  float *Pop_O;

  Pop_O = vector(1,n);
  m = 0;
  for (i=0; i<n; i++)
    if(Pop[i] > 0)
      Pop_O[++m] = -Pop[i];
  sort(m,Pop_O);
  /* revert(m+1,Pop_O); */
  for (i=0; i<m; i++) Rank[i] += -(double)Pop_O[i+1];
  //for (i=m; i<n; i++) Rank[i] += 0.;
  /* Rank[0] is storing the maximum value of the n-dimensional vector
     (Pop[0],...,Pop[n-1]), where Rank[i=n] would be actually always equal to zero,
     if Rank[i] were defined from i=0 to i=n.
  */
  free_vector(Pop_O, 1,n);
}

void ordering_Vector_double(double *X, int n, int HOW)
{
  /* if HOW = 1, the output vector is ordered in increasing values */
  /* if HOW =-1, the output vector is ordered in decreasing values */
  int i,m;
  float *Pop_O;

  Pop_O = vector(1,n);
  m = 0;
  for (i=0; i<n; i++)
    Pop_O[++m] = (float)HOW *(float)X[i];  /* Pop_O[++m] = (float)X[i]
                                                Output vector ordered in
						increasing values */
  sort(m,Pop_O);
  for (i=0; i<n; i++)
      X[i] = (double)HOW * (double)Pop_O[i+1]; /* X[i] = (double)Pop_O[i+1]
				                Output vector odered in
						increasing values */
  free_vector(Pop_O, 1,n);
}

void ordering_Vector_2_double(double *X, double *Y, int n, int HOW)
{
  /* Sorts arrays X into ascending order using Quicksort, while mantaining the
     corresponding rearrangement of the array Y */
  /* if HOW = 1, the output vector is ordered in increasing values */
  /* if HOW =-1, the output vector is ordered in decreasing values */
  int i,m;
  float *Pop_O, *Pop_A;

  Pop_O = vector(1,n); Pop_A = vector(1,n);
  m = 0;
  for (i=0; i<n; i++){
    m++;
    Pop_O[m] = (float)HOW *(float)X[i];  /* Pop_O[++m] = (float)X[i]
                                                Output vector ordered in
						increasing values */
    Pop_A[m] = (float)Y[i];
  }
  sort2(m,Pop_O,Pop_A);
  for (i=0; i<n; i++){
      X[i] = (double)HOW * (double)Pop_O[i+1]; /*
						  X[i] = (double)Pop_O[i+1]
				                  Output vector odered in
						  increasing values
					       */
      Y[i] = (double)Pop_A[i+1];
  }
  free_vector(Pop_O, 1,n); free_vector(Pop_A, 1,n);
}

void rank_Abundance_Ave(double Pop[], int n, double Rank[], double Rank_Var[])
{
  int i,m;
  float *Pop_O;
  double *Rank_new;

  Rank_new = dvector(0,n-1);
  Pop_O = vector(1,n);
  m = 0;
  for (i=0; i<n; i++)
    if(Pop[i] > 0.)
      Pop_O[++m] = -(float)Pop[i];
  sort(m,Pop_O);
  /* revert(n+1,Pop_O); */
  for (i=0; i<m; i++) Rank_new[i] = -(double)Pop_O[i+1];
  for (i=m; i<n; i++) Rank_new[i] = 0.;
  Norma_Total_0_J(Rank_new,n-1);

  for (i=0; i<n; i++){
    Rank[i] += Rank_new[i];
    Rank_Var[i] += Rank_new[i] * Rank_new[i];
  }
  free_dvector(Rank_new, 0,n-1);
  free_vector(Pop_O, 1,n);
}
#endif

void Ave_Var(double *Ave, double *Var, int Max_Sp, int n)
{
  int i;

  for(i=0; i<Max_Sp; i++){
    Ave[i] /= (double)n;
    Var[i] /= (double)n;
    Var[i] = sqrt(Var[i] - Ave[i]*Ave[i]);
  }
}

void accum_Histogram(double *Sim, int No, double *y, double *y_V)
{
  int i;

  for(i=0; i<No; i++){
    y[i] += Sim[i];
    y_V[i] += Sim[i]*Sim[i];
  }
}

void accum_Diversity_Histogram(int *Divers, int No, double *y, double *y_V)
{
  int i;

  for(i=1; i<=No; i++){
    y[i-1] += (double)Divers[i];
    y_V[i-1] += (double)Divers[i]*(double)Divers[i];
  }
}

void rank_Abundace_Ave_Var(double Rank_Ave[], double Rank_Var[], int Max_Sp, int n)
{
  int i;

  for(i=0; i<Max_Sp; i++){
    Rank_Ave[i] /= (double)n;
    Rank_Var[i] /= (double)n;
    Rank_Var[i] = sqrt(Rank_Var[i] - Rank_Ave[i]*Rank_Ave[i]);
  }
}

double vector_MAX_double(double * P, int n_0, int No)
{
  /* MAX value among P[n_0], ..., P[No-1] */
  int i;
  double F, max_Value;

  F = P[n_0];
  for(i=n_0; i<No; i++){
    F = MAX(P[i],F);
  }
  max_Value = F;

  return (max_Value);
}

double vector_MAX_Index_double(double * P, int n_0, int No, int * I)
{
  /* MAX value among P[n_0], ..., P[No-1] */
  int i;
  double F, max_Value;

  F = P[n_0];
  for(i=n_0; i<No; i++){

    if ( F <= P[i] ) {
      (*I) = i;
      F = P[i];
    }

  }

  max_Value = F;

  return (max_Value);
}

int vector_MAX_int(int * P, int n_0, int No)
{
  /* MAX value among P[n_0], ..., P[No-1] */
  int i;
  int F, max_Value;

  F = P[n_0];
  for(i=n_0; i<No; i++){
    F = MAX(P[i],F);
  }
  max_Value = F;

  return (max_Value);
}


double vector_MIN_double(double * P, int n_0, int No)
{
  /* MIN value among P[n_0], ..., P[No-1] */
  int i;
  double F, min_Value;

  F = P[n_0];
  for(i=n_0; i<No; i++){
    F = MIN(P[i],F);
  }
  min_Value = F;

  return (min_Value);
}

int discrete_MAX_double(double *P, int No)
{
  int i, i_max;
  double F;

  F = P[0];
  for(i=0; i<No; i++){
    if(P[i] >= F) i_max = i;
    F = MAX(P[i],F);
  }

  return i_max;
}

/* Functions to save a vector in a predefined, ordered file */
void Saving_to_File_int(char Prefix[], double Time[], int Number[],
		    int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
  */
  int i;
  char file[12];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\t %d\n",Time[i], Number[i]);
  fclose(fp);
}

void name_Ordered(char Prefix[], int File_Number, char Sufix[], char *File)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
  */
  int i;

  File[0]='\0';
  fitxer(File, Prefix, File_Number, Sufix);
  printf("Generated File name: %s\n",File);
}

void Saving_to_File_distribution(char Prefix[], double Ps[],
				 int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
  */
  int i;
  char file[40];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving Distribution in the file %s\n",file);
  fp = fopen(file, "w");
  for(i=0; i < Points; i++)
    if(Ps[i] > 0.0)
      fprintf(fp,"%d\t%g\n",i,Ps[i]);
  fclose(fp);
  printf("Success: Distribution has been saved in the file %s!\n",file);
}

void Saving_to_File_float_float(char Prefix[], float Time[], float Number[],
				int Points, int File_Number, int ZEROS)
{
  int i;
  char file[24];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s ... ...\t",file);
  fp = fopen(file, "w");
  for(i=0; i<Points; i++){
    if(ZEROS == 1){
      if(Number[i] > 0.)  fprintf(fp,"%f\t%f\n", Time[i], Number[i]);
    }
    else {
      fprintf(fp,"%f\t%f\n", Time[i], Number[i]);
    }
  }
  fclose(fp);

  printf("File %s has been successfully saved!!!\n",file);
}

void Saving_to_File_float(char Prefix[], double Time[], float Number[],
			  int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
  */
  int i;
  char file[12];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\t %f\n",Time[i], Number[i]);
  fclose(fp);
}

void saving_Time_Two_Float_vectors(char *file,
				  float *time, float *x, float *y, int nP)
{
  int i;
  FILE *fp;

  fp = fopen(file, "w");
  printf("Saving %s... \n", file);
  for(i = 0; i<nP; i++){

    fprintf(fp, "%f\t%f\t%f\n",time[i],x[i],y[i]);
  }
  fclose(fp);
}

void Reading_from_File_double_nxy_Counting_Rows(char * File_Name,
						double * x_Data, double ** y_Data,
						int MAX_No_of_ROWS,
						int n_XY,
						int * No)
{
  /* This function counts the nubmer of rows (int * No) of the file while to read
     It checks if the vector and Matrix to read will be too large for the quantity
     of memory reserved in main code (through MAX_No_of_ROWS)
  */
  int n, j;
  double x,y;
  FILE *fp;

  int Max_No_of_ROWS;

  printf("\n[Entering function Reading_from_File_double_nxy() from stat.c ]\n    Reading File %s...\n",
	 File_Name); //Press_Key();

  if((fp=fopen(File_Name,"r")) == NULL) {
    printf("File non-existent! Cannot open file.\n");
    exit(1);
  }

  n=0;
  while ( fscanf(fp, "%lf\t", &x) != EOF ){
    x_Data[n] = x;
    for(j=0; j < n_XY; j++){
      fscanf(fp, "%lf\t", &y);
      y_Data[j][n] = y;
    }

    n++;
    assert( n <= MAX_No_of_ROWS );
  }

  * No = n; /* Number of rows (for instance, times) in the file to read */

  fclose(fp);
  printf("    File %s has been read successfully\n", File_Name); //Press_Key();
  printf("[Exiting function Reading_from_File_double_nxy() from stat.c]\n\n");
}

void Reading_from_File_double_nxy(char * File_Name, double * x_Data, double ** y_Data, int No, int n_XY)
{
  int i,j;
  double x,y;
  FILE *fp;

  printf("\n[Entering function Reading_from_File_double_nxy() from stat.c ]\n    Reading File %s...\n",
	 File_Name); //Press_Key();
  if((fp=fopen(File_Name,"r")) == NULL) {
    printf("File non-existent! Cannot open file.\n");
    exit(1);
  }
  for(i=0; i < No; i++){
    fscanf(fp, "%lf\t", &x);
    x_Data[i] = x;
    for(j=0; j < n_XY; j++) { fscanf(fp, "%lf\t", &y); y_Data[j][i] = y; }
    fscanf(fp, "\n");
  }
  fclose(fp);
  printf("    File %s has been read successfully\n", File_Name); //Press_Key();
  printf("[Exiting function Reading_from_File_double_nxy() from stat.c]\n\n");
}

void Reading_from_File_double_y(char * File_Name, double ** y_Data, int No_Rows, int No_Cols)
{
  int i,j;
  double y;
  FILE *fp;

  printf("\n[Entering function Reading_from_File_double_y() from stat.c ]\n    Reading File %s...\n",
	 File_Name); //Press_Key();
  if((fp=fopen(File_Name,"r")) == NULL) {
    printf("File non-existent! Cannot open file.\n");
    exit(1);
  }
  for(i=0; i < No_Rows; i++){
    for(j=0; j < No_Cols; j++) {
      if ( j == (No_Cols-1) ) { fscanf(fp, "%lf\n", &y); y_Data[i][j] = y; }
      else                    { fscanf(fp, "%lf\t", &y); y_Data[i][j] = y; }
    }
    //for(j=0; j < No_Cols; j++){ printf("%g\t", y_Data[i][j]); }
    //Press_Key();
    //fscanf(fp, "\n");
  }
  fclose(fp);
  printf("    File %s has been read successfully\n", File_Name); //Press_Key();
  printf("[Exiting function Reading_from_File_double_y() from stat.c]\n\n");
}

void Saving_to_File_double_y(char * File_Name, double ** y_Data, int No_Rows, int No_Cols)
{
  int i,j;
  double y;
  FILE *fp;

  printf("\n[Entering function Saving_to_File_double_y() from stat.c ]\n  Saving File %s...\n",
	 File_Name); //Press_Key();

  fp=fopen(File_Name, "w");
  for(i=0; i < No_Rows; i++){
    for(j=0; j < No_Cols; j++) { fprintf(fp, "%g\t", y_Data[i][j]); }
    fprintf(fp, "\n");
  }

  fclose(fp);
  printf("    File %s has been saving successfully\n", File_Name); //Press_Key();
  printf("[Exiting function Saving_to_File_double_y() from stat.c]\n\n");
}

void Reading_from_File_double_TwoCol(char *file, double Time[], double Number[],
				     int Points)
{
  int i;
  FILE *fp;
  double x,y;

  printf("[Entering reading function]\n");
  printf(" Trying Reading File %s [ %d rows to read ] \n", file, Points);

  if((fp=fopen(file,"r")) == NULL) {
    printf("  File non-existent in current working directory. Cannot open file.\n");
    exit(1);
  }

  printf(" File %s [ %d rows] is there to read\n", file, Points);
  for(i=0; i<Points; i++){
    fscanf(fp,"%lf\t%lf\n",&x, &y);
    Time[i] = x;
    Number[i] = y;
    //////////////printf("%g\t%g\n",x, y);
  }
  fclose(fp);

  printf("[Exiting reading function in Reading_from_File_double_TwoCol() from stat.c]\n");

}

void Saving_to_File_Lattice_double(char *Prefix, double *x, double *y, double **Lattice, int Nx, int Ny,
				   int File_Number)
{
  /* The number of individuals inhabiting the lattice changes in number with time.
     This function saves an snabshot at a given time and creates a file according to
     File_Number-th label.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
  */
  int i,j;
  char file[40];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s\n",file);
  fp = fopen(file, "w");

  for(i=0; i<Nx; i++){
    for(j=0; j<Ny; j++){
      fprintf(fp,"%g\t%g\t%g\n", x[i], y[j], Lattice[i][j]);
    }
  }
  fclose(fp);

  printf("File %s saved successfully!\n",file);
}

void Saving_to_File_double(char *Prefix, double Time[], double Number[],
			   int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
     Two colomns saving file...
  */
  int i;
  char file[100];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s\n",file);
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\t%g\n",Time[i], Number[i]);
  fclose(fp);
}

void Saving_to_File_double_1COL(char *Prefix, double Number[], int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
     Two colomns saving file...
  */
  int i;
  char file[50];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s\n",file);
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\n", Number[i]);
  fclose(fp);
}

void Saving_to_File_double_1Col(char *Prefix, double Number[], int Points)
{
  /*
     One column file is saved...
  */
  int i;
  char file[40];
  FILE *fp;
  char *f;

  file[0]='\0';
  f = strcat(file, Prefix);
  f = strcat(file, ".dat");

  printf("Saving file %s\n",file);

  fp = fopen(file, "w"); for(i=0; i<Points; i++) fprintf(fp,"%g\n",Number[i]);  fclose(fp);
}

void Saving_to_File_double_Average_plusminus_VAR(char Prefix[], double Time[], double Number[],
						 double Var[], int Points, int File_Number)
{
  /* File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
     Three colomns saving file!!!!
  */
  int i;
  char file[50];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s\n",file);
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\t%g\t%g\t%g\n",Time[i], Number[i]-Var[i], Number[i], Number[i]+Var[i]);
  fclose(fp);
}

void Saving_to_File_double_3Col(char Prefix[], double Time[], double Number[],
				double Var[], int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
     Three colomns saving file!!!!
  */
  int i;
  char file[50];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s\n",file);
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\t%g\t%g\n",Time[i], Number[i], Var[i]);
  fclose(fp);
}

void Saving_to_File_double_4Col(char Prefix[],
				double Time[], double y_1[], double y_2[], double y_3[],
				int Points, int File_Number)
{
  /* The number of species inhabiting the lattice changes in number with time.
     This function saves the temporal evolution of this quantitity of the
     File_Number-th stochastic realization.
     Dependencies: Saving_to_File() -> fitxer() -> inttochar().
     Three colomns saving file!!!!
  */
  int i;
  char file[40];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, File_Number, ".dat");
  printf("Saving file %s\n",file);
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%g\t%g\t%g\t%g\n",Time[i], y_1[i], y_2[i], y_3[i]);
  fclose(fp);
}

void save_a_Matrix_to_File_int(char Prefix[],int No,int **InEvol, int NoRow,int NoCol)
{
  int i,j;
  char file[12];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, No, ".dat");

  fp = fopen(file, "w");
  for(i=0; i<NoRow; i++){
    fprintf(fp,"%d\t",i);
    for(j=0; j<NoCol; j++)
      fprintf(fp,"%d\t",InEvol[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void save_a_Matrix_to_File_float(char Prefix[],int No,float **InEvol, int NoRow,int NoCol)
{
  int i,j;
  char file[12];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, No, ".dat");

  fp = fopen(file, "w");
  for(i=0; i<NoRow; i++){
    fprintf(fp,"%d\t",i);
    for(j=0; j<NoCol; j++)
      fprintf(fp,"%f\t",InEvol[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void save_a_TimeMatrix_to_File_int(char Prefix[], int No, double Time[], int **InEvol,
				   int NoRow,int NoCol)
{
  int i,j;
  char file[12];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, No, ".dat");

  fp = fopen(file, "w");
  for(i=0; i<NoRow; i++){
    fprintf(fp,"%g\t",Time[i]);
    for(j=0; j<NoCol; j++)
      fprintf(fp,"%d\t",InEvol[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void save_a_TimeMatrix_to_File_float(char Prefix[],int No, double Time[], float **fEvol,
				     int NoRow,int NoCol)
{
  int i,j;
  char file[12];
  FILE *fp;

  file[0]='\0';
  fitxer(file, Prefix, No, ".dat");

  fp = fopen(file, "w");
  for(i=0; i<NoRow; i++){
    fprintf(fp,"%g\t",Time[i]);
    for(j=0; j<NoCol; j++)
      fprintf(fp,"%f\t",fEvol[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

/* Related functions on which Saving_to_file *** functions depend */
void fitxer (char file[], char a[], int b, char c[])
{
  /* This function returns the name of a file where 'char a[]' is the prefix,
     'int b' is an ordering number i 'char c[]' is a sufix.
     For instance, noSp_ +  5 + .dat = noSp_5.dat.
  */
  char *pFile;
  char num[12];

  pFile = strcat(file,a);
  inttochar(b,num);
  pFile = strcat(file,num);
  pFile = strcat(file,c);
}

void label_Name( char * Name_OUT, char * Name_IN )
{
  /* This functions assigns a name onto Name_OUT.
     Usually Name_IN is just "[name]" and Name_OUT
     nees to be allocated in memmory in advance
  */
  char * f;

  Name_OUT[0] ='\0';
  f = strcat( Name_OUT, Name_IN );
}

void fileName (char file[], char a[], char c[])
{
  /* This function returns the name of a file where 'char a[]' is the prefix,
     and 'char c[]' is a sufix. If a[] = 'Barcelona' and c[]='.dat'
     for instance, then barcelona + .dat = barcelona.dat would the output
  */
  char *pFile;

  pFile = strcat(file,a);
  pFile = strcat(file,c);
}

void inttochar(int b, char num[])
{
  int kount;
  int i, x, xif;
  int div;

  div=1; xif=1;
  kount = 0;
  while (xif != 0){
    div *= 10;
    xif = b/div;
    kount++;
  }                      /* kount, number of digits */
  for (i=kount;i>0;i--){
    div /= 10;
    x = b/div;
    b = b - x * div;
    num[kount-i] = x + '0';
  }
  num[kount] ='\0';
}

int doubletochar(double value, char * Number)
{
  /*
     . value is the double precision number to be
     converted

     . Number [] MUST be big enough to hold all
     the characters of your number!!

     . mm is the number of meaningful digits will retained
  */
  int i;


  i = sprintf(Number, "%lf", value);

  //i = sprintf(Number, "%f", value);

  return 0;
}

void Saving_floatVector(char* immVector, float *mVector, int Sp)
{
  int i;
  FILE *fp;

  fp = fopen(immVector, "w");

  for(i=0; i<Sp; i++)
    fprintf(fp, "%f\n", mVector[i]);

  fclose(fp);
}

void Saving_intVector(char* immVector, int *mVector, int Sp)
{
  int i;
  FILE *fp;

  fp = fopen(immVector, "w");

  for(i=0; i<Sp; i++)
    fprintf(fp, "%d\n", mVector[i]);

  fclose(fp);
}

void Reading_floatVector(char* File, float *mVector, int Sp)
{

  int i;
  FILE *pf;
  float D;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");
  mVector[0] = 0.;
  for(i=1; i<=Sp; i++){
    fscanf(pf, "%f\n", &D);
    mVector[i] = D;
  }
}

void Reading_doubleVector(char* File, double *mVector, int Sp)
{

  int i;
  FILE *pf;
  double D;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");

  if (!pf){
    printf(" File %s is not existing (o sea, brilla por su ausencia). The program will exit... ", File); Press_Key();
    exit(0);
    return;
  }

  mVector[0] = 0.;
  for(i=0; i<Sp; i++){
    fscanf(pf, "%lf\n", &D);
    mVector[i] = D;
  }
}

void Reading_intVector(char* File, int *mVector, int Sp)
{

  int i;
  FILE *pf;
  int D;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");

  if (!pf){
    printf(" File %s is not existing (o sea, brilla por su ausencia). The program will exit... ", File); Press_Key();
    exit(0);
    return;
  }

  mVector[0] = 0;
  for(i=0; i<Sp; i++){
    fscanf(pf, "%d\n", &D);
    mVector[i] = D;
  }
}

void Reading_float_iniVector(char* File, float *mVector, int Sp)
{

  int i;
  FILE *pf;
  float D;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");
  for(i=0; i<Sp; i++){
    fscanf(pf, "%f\n", &D);
    mVector[i] = D;
  }
}

int Sum_iVector(int *mVector, int Sp)
{

  int i, Summ;

  Summ = 0;
  for(i=0; i<Sp; i++)
    Summ += mVector[i];

  return Summ;
}

double Sum_dVector(double *mVector, int Sp)
{

  int i;
  double Summ;

  Summ = 0.;
  for(i=0; i<Sp; i++)
    Summ += mVector[i];

  return Summ;
}

int Sum_Diversity(int *mVector, int Sp)
{

  int i, Summ;

  Summ = 0;
  for(i=0; i<Sp; i++)
    Summ += mVector[i] * i;

  return Summ;
}

double Sum_Diversity_double(double *mVector, int Sp)
{

  int i;
  double Summ;

  Summ = 0;
  for(i=0; i<Sp; i++)
    Summ += mVector[i] * i;

  return Summ;
}

void Reading_DiversityData(char* File, int *Species, int Data)
{
  int i;
  FILE *pf;
  int n, Sp;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");
  for(i=0; i<Data; i++){
    fscanf(pf, "%d\t%d\n", &n, &Sp);
    Species[n] = Sp;
  }
  fclose(pf);
}

void Reading_DiversityData_double(char* File, int *Species, int Data)
{
  int i;
  FILE *pf;
  int n;
  double Sp;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");
  for(i=0; i<Data; i++){
    fscanf(pf, "%d\t%lf\n", &n, &Sp);
    Species[n] = (int)Sp;
  }
  fclose(pf);
}

void Reading_DiversityData_double_BIS(char* File, double *Species, int Data)
{
  int i;
  FILE *pf;
  int n;
  double Sp;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");
  for(i=0; i<Data; i++){
    fscanf(pf, "%d\t%lf\n", &n, &Sp);
    Species[n] = Sp;
  }
  fclose(pf);
}

void Reading_DiversityData_ONE_COLUMN(char* File, int *Species, int Data)
{
  int i;
  FILE *pf;
  int n, Sp;

  printf("Reading %s... \n", File);
  pf = fopen(File, "r");
  for(i=0; i<Data; i++){
    fscanf(pf, "%d\n", &Sp);
    Species[i] = Sp;
  }
  fclose(pf);
}
/* End of statistical functions */

void linearEqSystem_2dim(double *p0, double *p1, double *p2, double *Sol)
{
  /*
     The equation is given as:
     p0[0] = p1[0] * x + p2[0] * y;
     p0[1] = p1[1] * x + p2[1] * y;

     Here, x = Sol[0] and y = Sol[1];
  */

  Sol[0] = (p2[1]*p0[0] - p0[1]*p2[0])/(p1[0]*p2[1]-p1[1]*p2[0]);
  Sol[1] = (p1[1]*p0[0] - p0[1]*p1[0])/(p1[1]*p2[0]-p1[0]*p2[1]);
}

int Finding_a_Value_int(const int * N, int S, int value,  int * i_th)
{
  /* This function finds the i-th position where N[i] == value within
     a vector N of length S. It returns the number of times this
     value is found and the last position, i-th, where this actually
     happens.
  */
  int i;
  int n;

  n = 0;
  for (i=0; i<S; i++){

    if(N[i] == value){

      (* i_th) = i;
      n++;

    }
  }

  return(n);
}

void create_Title_MONTH_YEAR(char * Title, char * Pre_Title, int i, int Initial_Year, int Initial_Month)
{
  /* Enough memmory has been allocated for Title when Title is passed to this function */

  int j;
  char *pFile;
  char CURRENT_YEAR[5];
  char ** MONTH;
  char * CURRENT_MONTH;
  int Year, Month;

  MONTH = (char**)malloc( 12 * sizeof(char *) );
  for (j = 0; j<12; j++) {
    MONTH[j] = (char*)malloc( 3 * sizeof(char) );
    MONTH[j][0]='\0';
  }

  pFile = strcat(MONTH[0], "Jan");
  pFile = strcat(MONTH[1], "Feb");
  pFile = strcat(MONTH[2], "Mar");
  pFile = strcat(MONTH[3], "Abr");
  pFile = strcat(MONTH[4], "May");
  pFile = strcat(MONTH[5], "Jun");
  pFile = strcat(MONTH[6], "Jul");
  pFile = strcat(MONTH[7], "Aug");
  pFile = strcat(MONTH[8], "Set");
  pFile = strcat(MONTH[9], "Oct");
  pFile = strcat(MONTH[10], "Nov");
  pFile = strcat(MONTH[11], "Dec");

  for (j = 0; j<12; j++) {
    CURRENT_MONTH = MONTH[j];
    printf("%s ", CURRENT_MONTH);
  }
  printf("\n");

  Year  = Initial_Year  + i/12;
  Month = Initial_Month + i%12;

  inttochar(Year, CURRENT_YEAR);
  CURRENT_MONTH = MONTH[Month];

  Title[0]='\0';
  pFile = strcat(Title, Pre_Title);
  pFile = strcat(Title, CURRENT_YEAR);
  pFile = strcat(Title, " - ");
  pFile = strcat(Title, CURRENT_MONTH);

  for (j = 0; j<12; j++) free(MONTH[j]); free (MONTH);
}

int Number_of_Octaves( int max_A, double * x, int MAX_NUMBER_of_OCTAVES_EVER)
{
      int n_Bin;
      int i, j, Inf, Sup;

      Inf = 1; Sup = 1; i=0;
      while(Sup < max_A){
	Inf = Sup;
	Sup = 2 * Sup;

	if (i < MAX_NUMBER_of_OCTAVES_EVER)
	  x[i] = sqrt((double)Inf * (double)(Sup-1));
	else
	  {
	    printf(" MAX_NUMBER_of_OCTAVES_EVER = %d\ti = %d\n", MAX_NUMBER_of_OCTAVES_EVER, i);
	    printf(" Chose a larger MAX_NUMBER_of_OCTAVES_EVER from MODEL.h!!!\n");
	    printf(" You were about to go over range of array x[] in SAD_Control.c\n");
	    exit(0);
	  }
	i++;
      }

      n_Bin = i;

      return(n_Bin); /* Number of Bins in the plot representation */
}
