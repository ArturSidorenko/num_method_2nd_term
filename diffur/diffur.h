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
#include<algorithm>


//fun fills the array ans with f(x, y) for dimension N (N is needed for checks)
//if N is not of appropriate value, an std::invalid_value exeption throws
void fun(double *ans, const double *y, double x, int N);

void fexact(double *y, double x, int N);


//performs the step 
void step(double *ynew, const double *y, double x, double h, int N);

//solves the ode itself and allocates memory for the task
void solve(double ***ans, size_t *points, double h, int N, double x0, const double* y0, double xend);

//other versions to carry out solution via two previous values of y
void new_step(double *ynew, const double *y, const double *yold, double x, double h, int N);
void new_solve(double ***ans, size_t *points, double h, int N, double x0, const double* y0, double xend); 

//awesome writing
void write(const char *name, int N, int m, double x0, double h, double **y);

//Linf norm between the our solution and the exact one
double Linf_norm(int N, int m, double x0, double h, double **y);


