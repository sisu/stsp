#include "glptsp.hpp"
#include "cut.hpp"
#include <cstdlib>
#include <iostream>
using namespace std;

GLPTSP::GLPTSP()
{
	lp = 0;
	relaxationOnly = 1;
	first = 0;
}

GLPTSP::~GLPTSP()
{
	if (!lp) return;
	delete[] cols;
	delete[] row;
	delete[] vars;
	glp_delete_prob(lp);
}

void GLPTSP::init()
{
	int N = dists.size();
	int M = N*(N-1)/2;
	cols = new int[M+1];
	row = new double[M+1];
	vars = new double[M+2];
	lp = glp_create_prob();
	glp_add_cols(lp, M);
	for(int i=0; i<M; ++i) {
//		glp_set_col_kind(lp, 1+i, GLP_BV);
		if (i==0) glp_set_col_bnds(lp, 1+i, GLP_FX, 1, 1);
		else if (relaxationOnly) glp_set_col_bnds(lp, 1+i, GLP_DB, 0, 1);
		else glp_set_col_kind(lp, 1+i, GLP_BV);
	}

	for(int i=0; i<N; ++i) {
		int v = varnum(i, (i+1)%N);
		glp_set_col_stat(lp, v, GLP_BS);
	}

	glp_set_obj_dir(lp, GLP_MIN);
	int K=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			glp_set_obj_coef(lp, ++K, dists[i][j]);

	int r = glp_add_rows(lp, N);
	for(int i=1; i<=N; ++i) row[i] = 1;
	for(int i=0; i<N; ++i, ++r) {
		makeDegreeConstraint(cols+1, i, N);
		glp_set_mat_row(lp, r, N-1, cols, row);
		glp_set_row_bnds(lp, r, GLP_FX, 2, 2);
		glp_set_row_stat(lp, r, GLP_NL);
	}

	// TODO: better basis?
	glp_std_basis(lp);

	first = 1;
}
void GLPTSP::reset()
{
	int N = dists.size();
	int M = N*(N-1)/2;
	for(int i=0; i<M; ++i) {
//		glp_set_col_kind(lp, 1+i, GLP_BV);
		if (i==0) glp_set_col_bnds(lp, 1+i, GLP_FX, 1, 1);
		else if (relaxationOnly) glp_set_col_bnds(lp, 1+i, GLP_DB, 0, 1);
		else glp_set_col_kind(lp, 1+i, GLP_BV);
	}
}

static vector<int> arows;
double GLPTSP::calc()
{
	int N = dists.size();
	int M = N*(N-1)/2;
	int K=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			glp_set_obj_coef(lp, ++K, dists[i][j]);

	arows.clear();
	double r=0;
	if (relaxationOnly) r = relaxation(M);
	else r = branchAndCut(M);
	if (!arows.empty()) glp_del_rows(lp, arows.size(), &arows[0]-1);
	return r;
}

double GLPTSP::relaxation(int M)
{
	int N = dists.size();
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;
	parm.meth = GLP_PRIMAL;
	parm.presolve = first;
	first=0;

//	bool fst=1;
	do {
#if 0
		if (!fst) {
			cout<<"reoptimize\n";
		} else cout<<"fst "<<N<<'\n';
		fst=0;
#endif

		int r = glp_simplex(lp, &parm);
		if (r) {
			cout<<"FAILED SOLVING: "<<r<<'\n';
			abort();
		}
		for(int i=1; i<=M; ++i) vars[i] = glp_get_col_prim(lp, i);
		parm.meth = GLP_DUAL;
		parm.presolve = 0;

//		genSubtours(N, vars, *this);
//		break;
	} while(genSubtours(N, vars, *this));
	return glp_get_obj_val(lp);
}
double GLPTSP::branchAndCut(int M)
{
	glp_iocp parm;
	glp_init_iocp(&parm);
	parm.cb_func = GLPTSP::callback;
	parm.cb_info = this;
	parm.msg_lev = GLP_MSG_ERR;
	parm.tm_lim = 2500;
//	parm.presolve = GLP_ON;

#if 0
	glp_smcp pp;
	glp_init_smcp(&pp);
	pp.msg_lev = GLP_MSG_ERR;
	int rr = glp_simplex(lp, &pp);
	if (rr) {
		cout<<"SOLVING RELAXATION FAILED: "<<rr<<'\n';
		abort();
	}
#else
	double r0 = relaxation(M);
#endif
	int r = glp_intopt(lp, &parm);
	if (r==GLP_ETMLIM) {
		cout<<"time limit reached "<<glp_get_obj_val(lp)<<' '<<r0<<'\n';
		return glp_get_obj_val(lp);
	} else if (r) {
		cout<<"FAILED SOLVING INTOPT: "<<r<<'\n';
		abort();
	}
	cout<<"ratio: "<<r0 / glp_mip_obj_val(lp)<<'\n';
	return glp_mip_obj_val(lp);
}

void GLPTSP::addCut(int n, int* cols, double* row, double low)
{
	int r = glp_add_rows(lp, 1);
	glp_set_mat_row(lp, r, n, cols-1, row-1);
	glp_set_row_bnds(lp, r, GLP_LO, low, low);
	arows.push_back(r);
}

struct GLPWrap : LP {
	glp_prob* lp;
	GLPWrap(glp_prob* lp): lp(lp) {}
	void addCut(int n, int* cols, double* row, double low) {
		int r = glp_add_rows(lp, 1);
//		cout<<"adding row "<<r<<' '<<n<<' '<<glp_get_num_cols(lp)<<'\n';
		glp_set_mat_row(lp, r, n, cols-1, row-1);
		glp_set_row_bnds(lp, r, GLP_LO, low, low);
	}
};

void GLPTSP::callback(glp_tree* tree, void* info)
{
	switch(glp_ios_reason(tree)) {
		case GLP_IROWGEN: {
				GLPTSP* tsp = static_cast<GLPTSP*>(info);
				int N = tsp->dists.size();
				int M = N*(N-1)/2;
				glp_prob* lp = glp_ios_get_prob(tree);
//				cout<<"lol "<<N<<' '<<M<<' '<<glp_get_num_cols(lp)<<'\n';
				for(int i=1; i<=M; ++i) tsp->vars[i] = glp_get_col_prim(lp, i);
				GLPWrap w(lp);
				genSubtours(tsp->dists.size(), tsp->vars, w);
			}
			break;
		default:
			break;
	}
}
