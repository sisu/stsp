#include <vector>
#include <iostream>
#include "concorde.hpp"
#include "util.hpp"
#include "ctsp.hpp"
using namespace std;

extern vector<vector<double> > dist;

//typedef ConcordeTSP TSP;
typedef CustomTSP TSP;

namespace {

vector<TSP> tsps;

vector<vector<int> > samples;
vector<double> probs;

}

typedef vector<int> ivec;
double expectedTotalCost(const ivec& path)
{
	double r=0;
	for(size_t i=1; i<path.size(); ++i) {
		r += dist[path[i-1]][path[i]];
	}
	r *= .5;

	vector<double> pdist;
	for(size_t a=0; a<samples.size(); ++a) {
		vector<int>& v = samples[a];
		int n = v.size();
		pdist.resize(n);
		for(int i=0; i<n; ++i) {
			int a = v[i];
			double d=1e100;
			for(size_t j=0; j<path.size(); ++j)
				d = min(d, dist[a][path[j]]);
			pdist[i] = d;
		}

		for(int i=0; i<n; ++i)
			for(int j=0; j<i; ++j) {
				int x=v[i], y=v[j];
//				cout<<"asdasd "<<tsps[a].dists<<'\n';
				tsps[a].dists[i][j] = tsps[a].dists[j][i] = min(dist[x][y], pdist[i]+pdist[j]);
			}

		r += probs[a] * tsps[a].calc();
	}
	return r;
}

void initTSPCost(const vector<vector<int> >& ss, const vector<double> ps)
{
	samples.resize(ss.size());
	for(size_t i=0; i<samples.size(); ++i) {
		vector<int>& v = samples[i];
		v.resize(2+ss[i].size());
		v[0]=0;
		v[1]=1;
		for(size_t j=0; j<ss[i].size(); ++j)
			v[2+j] = 2+ss[i][j];
	}
	probs = ps;

	tsps.resize(ss.size());
	for(size_t i=0; i<tsps.size(); ++i) {
		TSP& tsp = tsps[i];
		int n = samples[i].size();
		cout<<"init tsp "<<i<<" : "<<n<<'\n';
		tsp.dists.resize(n, vector<double>(n));
		tsp.init();
	}
}
