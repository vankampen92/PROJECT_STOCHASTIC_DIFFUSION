#include <math.h>
#include <assert.h>
#include <nr.h>
#include <nrutil.h>
#define ITMAX 100
#define EPS 3.0e-7   /* Num Recipe EPS 3.0e-7 */
#define FPMIN 1.0e-30

double gammQ_double(double a, double x);
double gammP_double(double a, double x);
double gcf_D(double a, double x); 
double G_6_2_7(double a, double x);
double G_6_2_5(double a, double x); 
double gser_D(double a, double x); 
double Conditional_Gamma_Density_Distribution(double t, double n, double g); 

double gammQ_double(double a, double x)
{
	double gcf_D(double a, double x);
	double gser_D(double a, double x);
	void nrerror(char error_text[]);
	double gamser,gammcf;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gamser = gser_D(a,x);
		return (1.0-gamser);
	} else {
		gammcf = gcf_D(a,x);
		return gammcf;
	}
} 

double gammP_double(double a, double x)
{
	double gcf_D(double a, double x);
	double gser_D(double a, double x);
	void nrerror(char error_text[]);
	double gamser,gammcf;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gamser = gser_D(a,x);
		return gamser;
	} else {
		gammcf = gcf_D(a,x);
		return (1.0-gammcf);
	}
}

double gcf_D(double a, double x)
{
  double X;
  float gammln(float xx);
  float gln;
  double h;
  
  gln = gammln((float)a);
  
  h = G_6_2_7(a,x);
  
  X = exp(-x+a*log(x) - gln) * h;

  return(X);
}

double G_6_2_7(double a, double x)
{
  double h;
  void nrerror(char error_text[]);
  int i;
  double an,b,c,d,del;

  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");

  return (h);
}

double gser_D(double a, double x)
{
  double X;
  float gln;
  float gammln(float xx);

  gln = gammln((float)a);
  
  X = G_6_2_5(a, x) * exp(-x + a*log(x) - gln);
  
  return(X);
}

double G_6_2_5(double a, double x)
{
  void nrerror(char error_text[]);
  int n;
  double sum, sum_Fi, del, ap;

  if (a <= 0.0) nrerror("a less or equal to zero in routine G_6_2_5");
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    sum = 0.0;
    return (sum);
  } 
  else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	sum_Fi = sum;
	return (sum_Fi);
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return (sum);
  }
}
#undef ITMAX
#undef EPS
#undef FPMIN
	
double Conditional_Gamma_Density_Distribution(double t, double n, double g)
{
  /* Probability of recovery between time t and t+dt if the individual 
     was still infected at time t:
             Pr{ t < T < t+dt | Individual Still Infected at time t } 
     given that  T, the length of the infectious stage, is a stochastic variable
     distributed following a gamma distribution of (n,g) parameters, where
     1/g is the average length of the infectious stage.
  */
  double Ga,Z, x,X;
  double gln;
  float gammln(float xx);

  x = t*n*g;
  if(x < (n+1.0)){
    /* For times less than the mean infectious period */
    gln = gammln(n);
    X  = exp(n*log(x) - x - gln);
    Z  = X / (1. - X * G_6_2_5(n,x)); assert(Z>=0.0);
    Ga = 1./t * Z;
  }
  else {
    Z = -log(t * G_6_2_7(n,x));
    Ga = exp(Z);
  }
  
  return (Ga);
} 




