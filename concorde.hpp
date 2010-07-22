#ifndef CONCORDE_HPP
#define CONCORDE_HPP

#include <vector>

extern "C" {
#define new __new_tmp__
#include <tsp.h>
#undef new
}

struct ConcordeTSP {
	void init();
	double calc();
	std::vector<std::vector<double> > dists;

	bool relaxationOnly;

private:
	CCtsp_lp* tsp;
	CCtsp_lpcuts* pool;
	CCdatagroup* dat;
	CCrandstate* rstate;
	CCtsp_cutselect* sel;

	int* elist;
	int* elens;

	int* ptour;
};

#endif
