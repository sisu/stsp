#include "ctsp.hpp"
#include "util.hpp"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cassert>
using namespace std;

void CustomTSP::init()
{
}
typedef vector<int> ivec;

double CustomTSP::pathCost(const ivec& v) const
{
	double r=0;
	for(size_t i=1; i<v.size(); ++i)
		r += dists[v[i-1]][v[i]];
	return r;
}
static double randf()
{
	return rand()/(double)RAND_MAX;
}

double CustomTSP::calc()
{
	int N = dists.size();
	assert(N>=3);
	ivec cur(N);
	for(int i=0; i<N; ++i) cur[i]=i;
	swap(cur[1],cur[N-1]);

	ivec best = cur;
	double bdist = pathCost(cur);
	double cdist = bdist;
	for(double t=10; t>.1; t*=.9) {
		size_t a = 1+rand()%(N-3);
		size_t b = 1+rand()%(N-3);
		if (b>=a) ++b;
		else swap(a,b);

		int ca=cur[a], cb=cur[b];
//		cout<<"lol "<<a<<' '<<b<<' '<<ca<<' '<<cb<<'\n';
		double da = dists[cur[a-1]][cb] -  dists[cur[a-1]][ca];
		double db = dists[ca][cur[b+1]] - dists[cb][cur[b+1]];
		double d = da+db;
		if (d<0 || randf()<exp(-d/t)) {
			reverse(cur.begin()+a,cur.begin()+b+1);
			cdist = pathCost(cur);
			if (cdist < bdist) {
//				best = cur;
				bdist = cdist;
			}
		}
	}
	return bdist;
}
