#pragma once

#include<stdio.h>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<exception>
#include<algorithm>


/*
 * Eigenvalues and eigenvectors 
 * of equation y'' = - \lambda y.
 * 
 * The task is to check a proposed solution
 */


double eigenvalue(unsigned N, unsigned q);
double eigenvector(unsigned N, unsigned q, unsigned k);

double scal_product(unsigned N, unsigned p, unsigned q);
double left_hand(unsigned N, unsigned q, unsigned k);
double right_hand(unsigned N, unsigned q, unsigned k);

//residual for a particular equation including boundaries
double eq_res(unsigned N, unsigned q);

//residual for all equations and boundaries
double eq_res_all(unsigned N);


//checks difference between Gramian and identity matrices
double orthogonality_check(unsigned N);




