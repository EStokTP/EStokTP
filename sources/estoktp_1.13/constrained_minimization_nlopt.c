#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nlopt.h>

static int n;
//static const nlopt_algorithm algorithm_outer = NLOPT_AUGLAG_EQ;
static const nlopt_algorithm algorithm_outer = NLOPT_LD_SLSQP;

static const nlopt_algorithm algorithm_inner = NLOPT_LD_LBFGS;
//static const nlopt_algorithm algorithm_inner = NLOPT_LD_MMA;
//static const nlopt_algorithm algorithm_inner = NLOPT_LD_SLSQP;
//static const nlopt_algorithm algorithm_inner = NLOPT_G_MLSL_LDS;

static nlopt_result res;
static const double inner_tol = 1.0e-4;
static const double outer_tol = 1.0e-4;
//static const double inner_tol = 1.0e-5;
//static const double outer_tol = 1.0e-5;


#define F77RETTYPE int
#ifndef NOF77UNDERSCORE
#   define F77FUNNAME(fun) fun ## _
#else
#   define F77FUNNAME(fun) fun
#endif

static double *vbuf1, *vbuf2;


F77RETTYPE
F77FUNNAME(initial_guess) (double *q);

F77RETTYPE
F77FUNNAME(v1) (const double *q, double *res);

F77RETTYPE
F77FUNNAME(dv1) (const double *q, double *res);

F77RETTYPE
F77FUNNAME(v2) (const double *q, double *res);

F77RETTYPE
F77FUNNAME(dv2) (const double *q, double *res);

F77RETTYPE
F77FUNNAME(readv1) (double *res);

F77RETTYPE
F77FUNNAME(problem_size) (int *n);


static unsigned count = 0, gcount = 0, ccount = 0, gccount = 0;
double
myfunc (unsigned n, const double *x,
        double *grad, void *my_func_data)
{
  double V1val;

  F77FUNNAME(v1) (x, &V1val);
  ++count;  
  
  if (grad) {
    F77FUNNAME(dv1) (x, vbuf1);
    for (int ii = 0; ii < n; ++ii)
      grad[ii] = vbuf1[ii];
    gcount++;
  }
  //  for (int ii = 0; ii < n; ++ii)
  //    printf ("%d gradient V1 = %lf \n", ii, vbuf1[ii]);
  return V1val;
}

double
myconstfunc (unsigned n, const double *x,
             double *grad, void *my_func_data)
{
  double V1val;
  //  F77FUNNAME(v1) (x, &V1val);
  F77FUNNAME(readv1) (&V1val);

  double V2val;
  F77FUNNAME(v2) (x, &V2val);
  double retval;

  ccount++;
  retval = (V1val - V2val);
    
  if (grad) {
    F77FUNNAME(dv1) (x, vbuf1);
    F77FUNNAME(dv2) (x, vbuf2);
    for (int ii = 0; ii < n; ++ii)
      grad[ii] = (vbuf1[ii] - vbuf2[ii]);
    gccount++;
  }
  
  return retval;
}

//int na_tst_opt (void)
int na_tst_opt_()
{
  F77FUNNAME(problem_size) (&n);
  double *x = calloc (sizeof (double), n);
  vbuf1 = calloc (sizeof (double), n);
  vbuf2 = calloc (sizeof (double), n);
  printf ("problem size : %d\n", n);
  int nvar=n;
  FILE *nlopt_out;

  F77FUNNAME(initial_guess) (x);
  double f = HUGE_VAL;

  if((nlopt_out=fopen("nlopt_out.log","w"))==NULL) {
    printf("********failed opening %s*************\n","nlopt_out.log");
    exit(1);
  }

  printf ("first guess point %lf\n", x[0]);
  fprintf (nlopt_out,"first guess point %lf\n", x[0]);
  //  exit(0);
  
  nlopt_opt opt_outer = nlopt_create (algorithm_outer,  n);
  nlopt_opt opt_inner = nlopt_create (algorithm_inner,  n);

  res = nlopt_set_xtol_abs1 (opt_inner,  inner_tol);
  printf ("nlopt_set_xtol_abs completed with result %d\n", res);
  fprintf (nlopt_out,"nlopt_set_xtol_abs completed with result %d\n", res);

  //  res = nlopt_set_vector_storage (opt_inner, 100);
  //  printf ("nlopt_set_vector_storage completed with result %d\n", res);
  /*
  res = nlopt_set_local_optimizer (opt_outer, opt_inner);
  printf ("nlopt_set_local_optimizer completed with result %d\n", res);*/
  nlopt_destroy (opt_inner);
  
  res = nlopt_set_min_objective (opt_outer, myfunc, NULL);
  printf ("nlopt_set_min_objective completed with result %d\n", res);
  fprintf (nlopt_out,"nlopt_set_min_objective completed with result %d\n", res);

  res = nlopt_add_equality_constraint (opt_outer, myconstfunc, NULL, inner_tol);
  printf ("nlopt_add_equality_constraint completed with result %d\n", res);
  fprintf (nlopt_out,"nlopt_add_equality_constraint completed with result %d\n", res);

  res = nlopt_set_xtol_abs1 (opt_outer, outer_tol);
  printf ("nlopt_set_xtol_abs completed with result %d\n", res);
  fprintf (nlopt_out,"nlopt_set_xtol_abs completed with result %d\n", res);
  
  res = nlopt_optimize (opt_outer, x, &f);
  printf ("nlopt_optimize completed with result %d (a negative value means an error occurred)\n", res);
  fprintf (nlopt_out,"nlopt_optimize completed with result %d (a negative value means an error occurred)\n", res);

  if (res < 0)
    printf ("nlopt_optimize message : %s\n", nlopt_get_errmsg (opt_outer));
  fprintf (nlopt_out,"nlopt_optimize message : %s\n", nlopt_get_errmsg (opt_outer));
  
 
  printf ("optimum location : ");
  fprintf (nlopt_out,"optimum location : \n");
  for (int ii = 0; ii < n; ++ii)
    printf ("x[%d] = %g ", ii, x[ii]);
  for (int ii = 0; ii < n; ++ii)
    fprintf (nlopt_out,"x[%d] = %g \n", ii, x[ii]);

  printf ("\n");
  printf ("optimum value : f = %lf\n", f);
  printf ("number of function evaluations : %d\n", count);
  printf ("number of gradient evaluations : %d\n", gcount);
  printf ("number of constraint evaluations : %d\n", ccount);
  printf ("number of constraint gradient evaluations : %d\n", gccount);

  fprintf (nlopt_out,"\n");
  fprintf (nlopt_out,"optimum value : f = %g\n", f);
  fprintf (nlopt_out,"number of function evaluations : %d\n", count);
  fprintf (nlopt_out,"number of gradient evaluations : %d\n", gcount);
  fprintf (nlopt_out,"number of constraint evaluations : %d\n", ccount);
  fprintf (nlopt_out,"number of constraint gradient evaluations : %d\n", gccount);
  fclose(nlopt_out);
  nlopt_destroy (opt_outer);
  return 0;
}
