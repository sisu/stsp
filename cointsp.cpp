#include "cointsp.hpp"
#include "cut.hpp"
#include <iostream>
using namespace std;

CoinTSP::CoinTSP()
{
	cols=0;
	row=0;
	oneBased=0;
}
CoinTSP::~CoinTSP()
{
	delete[] cols;
	delete[] row;
}

void CoinTSP::init()
{
	lp.setLogLevel(0);
	int N = dists.size();
	int M = N*(N-1)/2;
	cols = new int[M+1];
	row = new double[M+1];
	lp.resize(0, M);
	for(int i=0; i<M; ++i) {
		if (i==0) lp.setColBounds(i, 1, 1);
		else lp.setColUpper(i, 1);
	}

	lp.setOptimizationDirection(1);
	int K=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			lp.setObjCoeff(K++, dists[i][j]);

	for(int i=0; i<N; ++i) row[i] = 1;
	for(int i=0; i<N; ++i) {
		makeDegreeConstraint(cols, i, N);
		lp.addRow(N-1, cols, row, 2, 2);
	}
}
void CoinTSP::reset()
{
}
double CoinTSP::calc()
{
	int N = dists.size();
	int K=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			lp.setObjCoeff(K++, dists[i][j]);

	arows.clear();
	double r=0;

	if (relaxationOnly) r = relaxation();
	else r = relaxation();

	if (!arows.empty()) lp.deleteRows(arows.size(), &arows[0]);
	return r;
}

double CoinTSP::relaxation()
{
	int N = dists.size();
	bool fst=1;
	double* vars=0;
	do {
		int r = fst ? lp.primal() : lp.dual();
		fst=0;
		if (r) {
			cout<<"FAILED SOLVING: "<<r<<'\n';
			abort();
		}
		vars = lp.primalColumnSolution();
	} while(genSubtours(N, vars, *this));
	return lp.getObjValue();
}

void CoinTSP::addCut(int n, int* cols, double* row, double low)
{
	int r = lp.getNumRows();
	arows.push_back(r);
	lp.addRow(n, cols, row, low);
}
