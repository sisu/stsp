#ifndef GLPTSP_HPP
#define GLPTSP_HPP

#include <vector>
#include <glpk.h>

struct GLPTSP {
	void init();
	double calc();
	std::vector<std::vector<double> > dists;

private:
	glp_prob* lp;
	int* cols;
	double* row;
};

#endif
