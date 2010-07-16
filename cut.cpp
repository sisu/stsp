#include <utility>
using namespace std;

int varnum(int a, int b)
{
	return a*(a-1)/2 + b + 1;
}

void makeDegreeConstraint(int* res, int k, int n)
{
	for(int i=0; i<k; ++i) res[i] = varnum(k, i);
	for(int i=k+1; i<n; ++i) res[i-1] = varnum(i, k);
}
