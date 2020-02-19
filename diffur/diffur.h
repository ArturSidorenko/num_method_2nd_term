/* Numerical resolving of a ODE y' = f(x, y)
 * via Euler method
 * where y and f are vectors
 */

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<exception>
#include<cstdlib>



//fun fills the array ans with f(x, y) for dimension N (N is needed for checks)
//if N is not of appropriate value, an std::invalid_value exeption throws
void fun(double *ans, double *y, double x, int N);

void fexact(double *y, double x, int N);


//performs the step 
void step(double *ynew, double *y, double x, double h, int N);

//solves the ode itself and allocates memory for the task
double** solve(double h, int N, double x0, double* y0, double xend);

//awesome writing
void write(const char *name, int N,  double x0, double xend, double h, double **y);

