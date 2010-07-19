#ifndef GLPTSP_HPP
#define GLPTSP_HPP

#include <vector>
#include <glpk.h>
#include "lp.hpp"

struct GLPTSP : LP {
	void init();
	double calc();
	std::vector<std::vector<double> > dists;
	void addCut(int n, int* cols, double* row, double low);

private:
	glp_prob* lp;
	int* cols;
	double* row;
	double* vars;
};

#endif
