#ifndef CUT_HPP
#define CUT_CPP

#include "lp.hpp"

int varnum(int a, int b);
void makeDegreeConstraint(int* res, int k, int n);
int genSubtours(int n, double* g, LP& lp);

#endif
