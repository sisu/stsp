#ifndef GLPTSP_HPP
#define GLPTSP_HPP

#include <vector>
#include <glpk.h>
#include "lp.hpp"

struct GLPTSP : LP {
	GLPTSP(): lp(0) {}
	~GLPTSP();
	void init();
	double calc();
	std::vector<std::vector<double> > dists;
	void addCut(int n, int* cols, double* row, double low);

private:
	glp_prob* lp;
	int* cols;
	double* row;
	double* vars;

	double relaxation(int M);
	double branchAndCut(int M);

	static void callback(glp_tree* tree, void* info);

	bool first;
};

#endif
