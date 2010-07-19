#include <utility>
#include <vector>
#include <iostream>
#include "lp.hpp"
using namespace std;

int varnum(int a, int b)
{
	if (a<b) swap(a,b);
	return a*(a-1)/2 + b + 1;
}

void makeDegreeConstraint(int* res, int k, int n)
{
	for(int i=0; i<k; ++i) res[i] = varnum(k, i);
	for(int i=k+1; i<n; ++i) res[i-1] = varnum(i, k);
}

namespace {
vector<double> conn;
vector<int> used;
vector<int> cols;
vector<double> row;
}

const double EPS = 1e-6;
int genSubtours(int n, double* g, LP& lp)
{
	conn.clear();
	conn.resize(n, 0);
	used.clear();
	used.resize(n, -1);

	row.clear();
	row.resize(n*(n-1)/2+1, 1);
	cols.reserve(row.size());

	int added = 0;
	for(int a=0; a<n; ++a) {
		if (used[a]>=0) continue;
		used[a]=a;
		for(int i=0; i<n; ++i) conn[i] = g[varnum(a,i)];
		while(1) {
			double bc=-1, bi=-1;
			for(int i=0; i<n; ++i) if (used[i]<0 && conn[i]>bc) bi=i, bc=conn[i];
			if (bc < .1) break;

			used[bi] = a;
			for(int i=0; i<n; ++i) conn[i] += g[varnum(bi,i)];

			double t=0;
			bool na=0;
			for(int i=0; i<n; ++i) if (used[i]!=a) t += conn[i], na=1;
			if (!na) break;
			if (t < 2 - EPS) {
				cols.clear();
				for(int i=0; i<n; ++i) if (used[i]==a)
					for(int j=0; j<n; ++j) if (used[j]!=a)
						cols.push_back(varnum(i,j));
				cout<<"addining cut "<<cols.size()<<' '<<t<<' '<<a<<' '<<n<<'\n';
				for(int i=0; i<n; ++i) cout<<used[i]<<' '; cout<<'\n';
				lp.addCut(cols.size(), &cols[0], &row[0], 2);
				break;
			}
		}
	}
//	cout<<"addeded "<<added<<" cuts\n";
	return added;
}
