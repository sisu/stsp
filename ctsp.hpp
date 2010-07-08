#ifndef CTSP_HPP
#define CTSP_HPP

#include <vector>

struct CustomTSP {
	void init();
	double calc();
	std::vector<std::vector<double> > dists;

private:
	double pathCost(const std::vector<int>& v) const;
};

#endif
