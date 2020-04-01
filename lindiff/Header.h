#pragma once

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<exception>
#include<cstdlib>
#include<algorithm>

#define MET1 1
#define MET2 2
#define MET3 3
#define MET4 4
#define MET5 5

void solve(double **ans, size_t *size, double A, double x0, double y0, double xe, double h, int nummet);

double step(double A, double y, double x, double h, int nummet);

double other_step(double A, double y, double y_old, double x, double h, int nummet);

double exact(double A, double x);
double L2_norm(double *y, size_t m, double h, double x0, double A);