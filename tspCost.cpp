#include <vector>
#include <iostream>
#include <algorithm>
#include "concorde.hpp"
#include "util.hpp"
#include "ctsp.hpp"
#include "glptsp.hpp"
#include "tspCost.hpp"
//#include "cointsp.hpp"
using namespace std;

extern vector<vector<double> > dist;
extern vector<vector<double> > edgeDist;
extern vector<vector<int> > conn;
extern bool robustOpt;

//typedef ConcordeTSP TSP;
//typedef CustomTSP TSP;
typedef GLPTSP TSP;
//typedef CoinTSP TSP;

namespace {

vector<TSP> tsps;

vector<vector<int> > samples;
vector<double> probs;

bool final=0;
}

typedef vector<int> ivec;
double expectedTotalCost(const ivec& path)
{
	double r=0;
	for(size_t i=1; i<path.size(); ++i) {
//		r += dist[path[i-1]][path[i]];
		int a = path[i-1], b = path[i];
		int n = lower_bound(conn[a].begin(),conn[a].end(),b)-conn[a].begin();
		assert(conn[a][n]==b);
		r += edgeDist[a][n];
//		if (final) cout<<"lol "<<a<<' '<<b<<' '<<n<<' '<<edgeDist[a][n]<<'\n';
	}
	if (final) cout<<path<<'\n'<<"length cost "<<r<<'\n';
	r *= LENGTH_FACTOR;

	double rr=0;

	vector<double> pdist;
	for(size_t k=0; k<samples.size(); ++k) {
		vector<int>& v = samples[k];
		int n = v.size();
		pdist.resize(n);
		for(int i=0; i<n; ++i) {
			int a = v[i];
			double d=1e100;
			for(size_t j=0; j<path.size(); ++j)
				d = min(d, dist[a][path[j]]);
			pdist[i] = d;
//			if (final) cout<<"pdist "<<k<<' '<<a<<": "<<d<<'\n';
		}

		for(int i=0; i<n; ++i)
			for(int j=0; j<i; ++j) {
				int x=v[i], y=v[j];
//				cout<<"asdasd "<<tsps[a].dists<<'\n';
				tsps[k].dists[i][j] = tsps[k].dists[j][i] = min(dist[x][y], pdist[i]+pdist[j]);
			}

		double cr = tsps[k].calc();
		if (robustOpt) rr = max(rr, cr);
		else rr += probs[k] * cr;

		if (final) cout<<"tsp cost "<<k<<": "<<cr<<'\n';
	}
	return r + rr;
}

double exactTotalCost(const ivec& path)
{
	final=1;
	for(size_t i=0; i<tsps.size(); ++i) {
		tsps[i].relaxationOnly = 0;
		tsps[i].reset();
	}
	return expectedTotalCost(path);
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
			v[2+j] = ss[i][j];
	}
	probs = ps;

	tsps.resize(ss.size());
	for(size_t i=0; i<tsps.size(); ++i) {
		TSP& tsp = tsps[i];
		int n = samples[i].size();
		cout<<"init tsp "<<i<<" : "<<n<<'\n';
//		for(int j=0; j<n; ++j) cout<<samples[i][j]<<' ';cout<<'\n';
		tsp.dists.resize(n, vector<double>(n));
		for(int j=0; j<n; ++j)
			for(int k=0; k<j; ++k)
				tsp.dists[j][k] = tsp.dists[k][j] = dist[samples[i][j]][samples[i][k]];
		tsp.init();
	}
}
