#ifndef COINTSP_HPP
#define COINTSP_HPP

#include "ClpSimplex.hpp"
#include "lp.hpp"
#include <vector>

struct CoinTSP : LP {
	CoinTSP();
	~CoinTSP();
	void init();
	void reset();
	double calc();
	std::vector<std::vector<double> > dists;
	void addCut(int n, int* cols, double* row, double low);
	bool relaxationOnly;

private:
	ClpSimplex lp;

	double relaxation();
	int* cols;
	double* row;
	std::vector<int> arows;
};

#endif
