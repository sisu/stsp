#ifndef LP_HPP
#define LP_HPP

struct LP {
	virtual void addCut(int n, int* cols, double* row, double low) = 0;
};

#endif
