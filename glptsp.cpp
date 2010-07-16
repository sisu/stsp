#include "glptsp.hpp"
#include "cut.hpp"
#include <cstdlib>
#include <iostream>
using namespace std;

void GLPTSP::init()
{
	int N = dists.size();
	int M = N*(N-1)/2;
	cols = new int[M+1];
	row = new double[M+1];
	lp = glp_create_prob();
	glp_add_cols(lp, M);
	for(int i=0; i<M; ++i) {
//		glp_set_col_kind(lp, 1+i, GLP_BV);
		if (i==0) glp_set_col_bnds(lp, 1+i, GLP_FX, 1, 1);
		else glp_set_col_bnds(lp, 1+i, GLP_DB, 0, 1);
	}

	glp_set_obj_dir(lp, GLP_MIN);
	int K=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			glp_set_obj_coef(lp, ++K, dists[i][j]);

	int r = glp_add_rows(lp, N+1);
	for(int i=1; i<=N; ++i) row[i] = 1;
	for(int i=0; i<N; ++i, ++r) {
		makeDegreeConstraint(cols+1, i, N);
		glp_set_mat_row(lp, r, N-1, cols, row);
		glp_set_row_bnds(lp, r, GLP_FX, 2, 2);
	}
	glp_std_basis(lp);
}

double GLPTSP::calc()
{
	int N = dists.size();
//	int M = N*(N-1)/2;
	int K=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			glp_set_obj_coef(lp, ++K, dists[i][j]);

	glp_smcp parm;
	glp_init_smcp(&parm);
	int r = glp_simplex(lp, &parm);
	if (r) {
		cout<<"FAILED SOLVING: "<<r<<'\n';
		abort();
	}
	return glp_get_obj_val(lp);
}
