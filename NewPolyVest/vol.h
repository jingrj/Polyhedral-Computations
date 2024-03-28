/***********************************************************************
*  This code is part of NewPolyVest.
***********************************************************************/

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "armadillo"
#include "time.h"
#include "memory.h"


#ifndef POLYVOL_H
#define POLYVOL_H

using namespace std;
using namespace arma;

namespace vol{

class polytope{
public:
	polytope(int rows, int cols);
	polytope(polytope& stu) {};
	polytope(){};
	~polytope();
	double 	matA(double val, int i, int j){ A(i, j) = val; return 0;}
	double 	matA(int i, int j){ return A(i, j); return 0;}
	double 	vecb(double val, int i){	b(i) = val; return 0;}
	double 	vecb(int i){ return b(i); return 0;}
	void 	Preprocess();
	double 	EstimateVol(int coef);
	void    Centering();
	rowvec  gaussrand();
	double  NewEstimateVol(int coef, double theta);
	
	
	int m, n;
	double determinant;
	
	bool 	msg_off;
	bool 	check_planes_off;
private:
	double 	walk(int k);
	void	checkHPs();
	void 	genInitE(double &R2, vec &Ori);

	double 	randd(double u){ return rand() * u / RAND_MAX; }
	int 	randi(int u){ return rand() % u; }

//polytope denoted by: Ax<=b, A is an (m x n) matrix.
	
	mat A;
	vec b;

//approximating volume variables
	vec x;
	double theta;
	FILE* data;
//	double vol;
	int l;
	double *r2;

	//reciprocal of beta, beta-cut
	double beta_r;

	vec *B;
	mat *Ai;
};

inline polytope::polytope(int rows, int cols) :
	msg_off(false),
	check_planes_off(false),
	m(rows),
	n(cols),
	A(rows, cols),
	b(rows),
	x(cols),
//	vol(0),
	determinant(0)
{
	srand((int)time(0));

	beta_r = 2 * n; //2 * n;

	l = (int)(n * log((double)beta_r) / log((double)2)) + 2;
	r2 = new double[l];
	for (int i = 0; i < l; i++) 
		r2[i] = pow((double)2, (double)(2 * i) / n);

	B = new vec[n];
	Ai = new mat[n];
}

inline polytope::~polytope(){
	delete []r2;
	delete []B;
	delete []Ai;
}
	
}

#endif

