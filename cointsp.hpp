#ifndef COINTSP_HPP
#define COINTSP_HPP

#include "ClpSimplex.hpp"

struct CoinTSP : LP {
	void init();
	void reset();
	double calc();
	std::vector<std::vector<double> > dists;
	void addCut(int n, int* cols, double* row, double low);
	bool relaxationOnly;

private:
	ClpSimplex lp;
};

#endif
