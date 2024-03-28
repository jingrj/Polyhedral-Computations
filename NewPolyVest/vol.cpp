/***********************************************************************
*  This code is part of NewPolyVest.
***********************************************************************/

#include <stdio.h>
#include <functional>
#include "vol.h"
#include "glpk.h"
#include <time.h>

#define PI 3.1415926536

using namespace vol;
using namespace std;
double uballVol(int n){
	double vol = 1;
	if (n % 2 == 1){
		int k = (n - 1) / 2;
		vol *= pow(2, n);
		for (int i = 1; i < k + 1; i++) vol *= PI * i;
		for (int i = 1; i < n + 1; i++) vol /= i;
	}else{
		int k = n / 2;
		for (int i = 1; i < k + 1; i++) vol *= PI / i;
	}
	return vol;
}

/*********** Delete Redundent Hyperplanes ***********/
void polytope::checkHPs(){
	if (check_planes_off) return;

	//init GLPK
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, m);
	glp_add_cols(lp, n);

	//disable msg output
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;

	//load constraints
	int *ind = new int[n + 1];
	double *val = new double[n + 1];
	for (int i = 1;i < m + 1; i++){
		for (int j = 1; j < n + 1; j++){
			ind[j] = j, val[j] = A(i - 1, j - 1);
		}
		glp_set_mat_row(lp, i, n, ind, val);
		glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
	}
	delete []ind, delete []val;
	for (int i = 1; i < n + 1; i++)
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

	//feasiblity check
	int num[2];
	for (int i = 1; i < m + 1;){
		glp_set_row_bnds(lp, i, GLP_LO, b(i - 1) + 0.00001, 0);
		glp_set_obj_coef(lp, 1, 1);
		for (int j = 1; j < n; j++)
			glp_set_obj_coef(lp, j + 1, 0);
		glp_simplex(lp, &parm);
		if (glp_get_status(lp) == GLP_NOFEAS){
			num[1] = i;
			glp_del_rows(lp, 1, num);
			A.shed_row(i - 1);
			b.shed_row(i - 1);
			m--;
		}else{
			glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
			i++;
		}
	}
	cout << "Hyperplanes Left: " << m << endl;
	glp_delete_prob(lp);
}

void polytope::genInitE(double &R2, vec &Ori){
	R2 = 0, Ori.zeros();

	//init GLPK
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, m);
	glp_add_cols(lp, n);

	//disable msg output
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;


	//load constraints
	int *ind = new int[n + 1];
	double *val = new double[n + 1];
	for (int i = 1; i < m + 1; i++){
		for (int j = 1; j < n + 1; j++){
			ind[j] = j, val[j] = A(i - 1, j - 1);
		}
		glp_set_mat_row(lp, i, n, ind, val);
		glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
	}
	delete []ind, delete []val;
	for (int i = 1; i < n + 1; i++)
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

	//get bounds
	for (int i = 0; i < n; i++){
		double max, min;
		for (int j = 0; j < n; j++)
			glp_set_obj_coef(lp, j + 1, 0);

		glp_set_obj_coef(lp, i + 1, 1);
		glp_simplex(lp, &parm);
		max = glp_get_obj_val(lp);
		for (int j = 0; j < n; j++)
			Ori(j) += glp_get_col_prim(lp, j + 1);

		glp_set_obj_coef(lp, i + 1, -1);
		glp_simplex(lp, &parm);
		min = -glp_get_obj_val(lp);
		for (int j = 0; j < n; j++)
			Ori(j) += glp_get_col_prim(lp, j + 1);

		R2 += (max - min) * (max - min);
	}
	Ori = Ori / (2 * n);
	
	glp_delete_prob(lp);
}

void polytope::Preprocess(){

	checkHPs();

	double c1 = (2 * pow(n, 2) + pow(1 - n / beta_r, 2)) * (1 - 1.0 / pow(beta_r, 2)) / (2 * pow(n, 2) - 2);
	//double c1 = pow(n, 2) * (1 - 1.0 / pow(beta_r, 2)) / (pow(n, 2) - 1);
	double c2 = (1 - n / beta_r) / (n + 1);
	double c3 = beta_r * beta_r;
	double c4 = 2 * c2 / (1 - 1.0 / beta_r);

	//init E(R2I, 0), T = R2I, ori = 0.
	mat T;
	vec ori(n);
	double R2;
	genInitE(R2, ori);
	T.eye(n, n);
	T = R2 * T;

	vec distance = zeros<vec>(m);
	vec tm = zeros<vec>(m);

	int counter = 0;
	while (++counter > 0){
		int i;
		
		//check if ori in polytope
		distance = b - A * ori;
		for (i = 0; i < m; i++)
			if (distance(i) < 0){
				tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
				break;
			}
		if (i == m){
			//check if small ellipsoid contained in polytope
			for (i = 0; i < m; i++){
				tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
				if (c3 * distance(i) * distance(i) - tm(i) < 0) break;
			}
		}
		
		//terminate if E satisfies two criterions
		if (i == m) break;
		
		vec t = T * A.row(i).t() / sqrt(tm(i));
		ori = ori - t * c2;
		T = c1 * (T - c4 * t * t.t());
	}
	
	if (!msg_off){ 
		//cout << "R^2: " << R2 << endl << "Origin: " << endl;
		//ori.print();
	}
	
	//apply affine transformation
	//mat Trans = chol(T);
	
	mat Trans;
	try{
		Trans = chol(T);
	}catch (const std::runtime_error& ex){
		cout << "The input polytope is degenerated or non-existed and the volume is 0." << endl;
		exit(1);		
	}
	
	//cout << Trans << endl;
	b = beta_r * (b - A * ori);
	A = A * Trans.t();

	//if (!msg_off) cout << "The number of iterations: " << counter << endl;

	rowvec exp(n);
	exp.ones();
	for (int i = 0; i < n; i++){
		B[i] = b / A.col(i);
		Ai[i] = A / (A.col(i) * exp);
	}
	
	determinant = det(Trans) / pow(beta_r, n);
}

void polytope::Centering() /* Find a better center point */
{
	vec newb = b;
	vec force = zeros(n,1);
	for (int j = 0; j < m; j++)
	{
		force -= norm(A.row(j)) / pow(b[j], 2) * A.row(j).t();
	}
	double t = 1;
	int count = 1;
	while (max(abs(force)) >= pow(10, -8) && count <= 1000)
	{
		newb = b - t * A * force; /* translate coordinate syste */
		if (min(newb) > 0) /* make sure the origin of the new coordinates system is still inside the polytope*/
		{
			vec newforce = zeros(n, 1);
			for (int j = 0; j < m; j++)
			{
				newforce -= norm(A.row(j)) / pow(newb[j], 2) * A.row(j).t();
			}
			if (norm(force) > norm(newforce))
			{
				b = newb; force = newforce; count++;
				continue;
			}
		}
		t = t / 2;
	}
}

rowvec polytope::gaussrand()
{
	rowvec allheight = zeros(1, m);
	for (int j = 0; j < m; j++)
	{
		allheight(j) = b(j)/norm(A.row(j));
	}
	arma::arma_rng::set_seed_random();
	vec r = randn(n, 1);
	r = r / norm(r);
	double distance = INFINITY;
	int rownum = 0;
	for (int j = 0; j < m; j++)
	{
		mat s = A.row(j) * r;
		if (s[0] > 0)
		{
			if (b[j] / s[0] < distance)
			{
				distance = b[j] / s[0];
				rownum = j;
			}
		}
	}
	rowvec sampledata = zeros(1, 2);
	sampledata[0] = distance;
	sampledata[1] = allheight[rownum];

	return sampledata;
}

