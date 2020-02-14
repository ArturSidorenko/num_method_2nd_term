/* Numerical resolving of a ODE y' = f(x, y)
 * via Euler method
 * where y and f are vector-valued fnctions
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

double exact_sol(double x, int entry);

//performs the step 
void step(double *ynew, double *y, double x, double h, int N);

//solves the ode itself and allocates memory for the task
double** solve(double h, int N, double x0, double y0, double xend);

//awesome writing
void write(const char *name, int N, int m, double x0, double h, double **y);

