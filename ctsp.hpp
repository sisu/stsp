#ifndef CTSP_HPP
#define CTSP_HPP

#include <vector>

struct CustomTSP {
	void init();
	double calc();
	void reset(){}
	std::vector<std::vector<double> > dists;
	bool relaxationOnly;

private:
	double pathCost(const std::vector<int>& v) const;
};

#endif
