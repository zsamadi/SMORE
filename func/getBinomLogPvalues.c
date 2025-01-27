/*=================================================================
 *
 * getBinomLogPvalues.C
 *
 * The calling syntax is:
 *
 *		[log_pvalue] = getBinomLogPvalues(t, y)
 *
 *
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"


/* Output Arguments */

#define log_pvalue_OUT plhs[0]


static double betacf(double a, double b, double x) {
 int MAXIT = 1000;
 double EPS = 3.0e-7;
 double FPMIN = 1.0e-300;
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  // These q's will be used in factors that occur in the coefficients
  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  // First step of Lentz's method
  c = 1.0;
  d = 1.0 - qab*x/qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1.0/d;
  h = d;
  for (m = 1; m <= MAXIT; m++) {
    m2 = 2*m;
    // One step (the even one) of the recurrence.
    aa = m*(b - m)*x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d*c;
    // Next step of the recurrence (the odd one).
    aa = -(a+m)*(qab + m)*x/((a + m2) * (qap + m2));
    d = 1.0 + aa*d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d*c;
    h *= del;
    if (fabs(del - 1.0) < EPS) break; // Are we done?
  }
  if (m > MAXIT) {
    fprintf(stderr, "a (%g) or b (%g) too big, or MAXIT (%d) too small in"
        " betacf. Input x is %g. Returning %g.\n",
        a, b, MAXIT, x, h);
  }
  return h;
}


static void getBinomLogPvalues(double log_pvalue[], double pos_count[], double neg_count[],double bernoulli, int m) {
		
	// (void)t; /* unused parameter */
	double log_pvaluet;

    double log_bt;
    int i;


	
	for (i=0; i<m; i++) {
		
        if (pos_count[i]>0.0){
    
          if (bernoulli < 0.0 || bernoulli > 1.0)  {
            fprintf(stderr, "Bad bernoulli in routine log_betai.\n");
            exit(EXIT_FAILURE);
          }
          if (bernoulli == 0.0 || bernoulli == 1.0) {
            log_bt = 1.0;
          } else {
            log_bt = lgamma(pos_count[i] +neg_count[i] + 1) - lgamma(pos_count[i]) - lgamma(neg_count[i] + 1)+ pos_count[i]*log(bernoulli)+ (neg_count[i]+1) *log(1.0 - bernoulli);
          }
          if (bernoulli < ((pos_count[i] + 1.0) / (pos_count[i] + (neg_count[i]+1) + 2.0))) {
            log_pvaluet=log_bt + log(betacf(pos_count[i], (neg_count[i]+1), bernoulli)/pos_count[i]);
          } else {
        
            log_pvaluet=log(1.0 - exp(log_bt) * betacf((neg_count[i]+1), pos_count[i], 1.0 - bernoulli) / (neg_count[i]+1));
          }
       

		  
		  log_pvalue[i]=log_pvaluet;
        }
        else {
            log_pvalue[i]=0;
        }


	
	}
    return;
}
 
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{

    double* log_pvalue;
    double  *pos_count, *neg_count, *bernoulli;


    size_t m;


    m = mxGetM(prhs[0]);


    /* Create a matrix for the return argument */
    log_pvalue_OUT = mxCreateDoubleMatrix((mwSize)m, (mwSize)1, mxREAL);

    /* Assign pointers to the various parameters */

    log_pvalue = mxGetPr(log_pvalue_OUT);
	


    pos_count = mxGetPr(prhs[0]);
    neg_count = mxGetPr(prhs[1]);
    bernoulli = mxGetPr(prhs[2]);
	
    /* Do the actual computations in a subroutine */
    getBinomLogPvalues(log_pvalue, pos_count,neg_count,bernoulli[0],m);
    return;
}

/* LocalWords:  log_pvalue maxlhs
 */

