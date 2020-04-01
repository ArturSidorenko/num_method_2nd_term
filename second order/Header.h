#pragma once

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<exception>
#include<cstdlib>
#include<algorithm>
#include<vector>


typedef double(*myfunc)(double);
typedef std::vector<double> myvec;


struct task {
	double h_;
	bool isconst_;
	double bconst_;
	myvec b_; 
	myvec f_;
	size_t size_;
	size_t N_;

	task(unsigned N, double b, myfunc right_hand);
	task(unsigned N, myfunc b, myfunc right_hand);
	double getb() const;
	double getb(size_t i) const;
};


//struct of a linear system
//the dimension equals m+1
//a, b, c encode a tridiagonal matrix
//c is the diagonal, a is the subdiagonal, b is the upperdiagonal
//indices of c vary from 0 to m
//indices of a vary from 1 to m
//indices of b vary from 0 to m-1
//indices of f vary from 0 to m
struct lin_sys {
	size_t m_;
	myvec a_, b_, c_;
	myvec f_;

	lin_sys(size_t m);
};

//linear system solvers
//boundary points are NOT included

//c is the diagonal, a is the subdiagonal, b is the upperdiagonal
//indices of c vary from 0 to m
//indices of a vary from 1 to m
//indices of b vary from 0 to m-1
//indices of f vary from 0 to m
//the answer is written into a vector v
void progonka(myvec & v, const lin_sys &s);
void fourier(myvec &v, const task &t);

//translate differential taks into a system of equations
void prepare_progonka(myvec &v, const task &t);

double basis(unsigned N, unsigned k, unsigned p); //a basis function
double eigenval(unsigned n, unsigned k);

double L2_norm(const myvec &a, const myvec &b, size_t m, double h);

