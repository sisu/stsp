#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "Vector.hpp"
#include "util.hpp"
using namespace std;

extern vector<vector<int> > conn;
extern int startI,endI;
extern vector<int> bestPath;
extern bool singleDir;

namespace {

vector<vector<double> > pheromone;
vector<vector<double> > probab;
vector<char> used;

double randf()
{
//	return fabs(rng()) / rng.max();
	return rand()/(double)RAND_MAX;
}

bool dfs(int s, int t, vector<int>& out)
{
	if (s==t) {
		out.push_back(s);
		return 1;
	}
//	cout<<"writing to "<<s<<'\n';
	int limit = 2-singleDir;
	++used[s];
	double p=0;
	for(size_t i=0; i<conn[s].size(); ++i) {
		int t = conn[s][i];
		if (used[t]<limit)
			p += probab[s][i];
	}
	double r = randf() * p;
	for(size_t i=0; i<conn[s].size(); ++i) {
		int x = conn[s][i];
		if (used[x]>=limit) continue;
		r -= probab[s][i];
		if (r>0) continue;
		if (dfs(x, t, out)) {
			out.push_back(s);
			return 1;
		} else {
			return dfs(s,t,out);
		}
	}
	return 0;
}

void addPheromone(const vector<int>& path, double p)
{
	for(size_t i=0; i+1<path.size(); ++i) {
		int a = path[i];
//		int e = find(conn[a].begin(),conn[a].end(),path[i+1]) - conn[a].begin();
		int e = lower_bound(conn[a].begin(),conn[a].end(),path[i+1]) - conn[a].begin();
		assert(e < (int)conn[a].size());
		pheromone[a][e] += p;
	}
}

}

double antColony(double(*cost)(const vector<int>&), double maxtime)
{
	const int T = 1<<11;
//	const int T = 1<<9;
	const int M = 1<<5;
	const double P = 5e-2;

	used.resize(conn.size());
	pheromone.resize(conn.size());
	probab.resize(conn.size());
	for(size_t i=0; i<conn.size(); ++i) {
		pheromone[i].resize(conn[i].size(),1);
		probab[i].resize(conn[i].size(),1);
	}
	double best = 1e100;

	int prevB = 0;
	vector<int> tmpv[M];
	double costs[M];

	clock_t start = clock();

	int a=0;
	for(a=0; a<T && a-prevB<T/4 && double(clock()-start)/CLOCKS_PER_SEC < maxtime; ++a) {
		int bnum=0;
		double bcost=1e100;
		for(int i=0; i<M; ++i) {
			fill(used.begin(),used.end(),0);
			vector<int>& path = tmpv[i];
			path.clear();
			dfs(startI, endI, path);
			reverse(path.begin(),path.end());

//			cout<<"lol path "<<path<<'\n';

			double c = costs[i] = cost(path);
//			cout<<"c "<<c<<'\n';
			if (c < bcost) {
				bcost=c;
				bnum = i;
			}
		}
		for(size_t i=0; i<pheromone.size(); ++i) {
			transform(pheromone[i].begin(),pheromone[i].end(),pheromone[i].begin(),bind1st(multiplies<double>(),1-P));
		}
		if (bcost<best) {
			best=bcost;
			bestPath = tmpv[bnum];
			cout<<"new best "<<best<<' '<<bestPath<<'\n';
			prevB=a;
		} else {
//			addPheromone(bestPath, .25*P);
			addPheromone(bestPath, 1./best);
		}
#if 0
		addPheromone(tmpv[bnum], .25*P);
#else
		for(int i=0; i<M; ++i) {
			addPheromone(tmpv[i], 1./costs[i]);
		}
#endif

		for(size_t i=0; i<probab.size(); ++i) {
			copy(pheromone[i].begin(),pheromone[i].end(),probab[i].begin());
		}
	}
	cout<<"antcolony ran for "<<a<<" iterations\n";

	return best;
}
