#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <cassert>
#include "util.hpp"
using namespace std;

namespace {
typedef pair<int,int> P;
typedef vector<int> ivec;

vector<ivec> nums;
vector<double> caps;

vector<double> exceeds;
vector<int> heights;
vector<int> destNode;

int N;

const double EPS = 1e-6;

void relabel(int n)
{
	int h = 1e9;
	for(size_t i=0; i<nums[n].size(); ++i) {
		int e = nums[n][i];
		if (caps[e]>EPS) h = min(h, heights[destNode[e]]);
	}
	heights[n] = 1+h;
}

void discharge(int n)
{
	size_t i=0;
	while(exceeds[n]>EPS) {
		if (i==nums[n].size()) {
			relabel(n);
			i=0;
		} else {
			int e = nums[n][i];
			int m = destNode[e];
//			if (m==0) cout<<"trying push: "<<heights[n]<<' '<<heights[m]<<' '<<caps[e]<<'\n';
			if (caps[e]>EPS && heights[m] < heights[n]) {
				double a = min(exceeds[n], caps[e]);
//				cout<<"pushing "<<n<<" -> "<<m<<" : "<<a<<'\n';
				exceeds[n] -= a;
				exceeds[m] += a;
				caps[e] -= a;
				caps[e^1] += a;
			}
			++i;
		}
	}
}

vector<char> used;
void resDFS(int n, vector<int>& res)
{
	used[n]=1;
	res.push_back(n);
	for(size_t i=0; i<nums[n].size(); ++i) {
		int e = nums[n][i];
		int t = destNode[e];
		if (caps[e]<EPS || used[t]) continue;
		resDFS(t, res);
	}
}

} // end anonymous namespace

double minCutFromTo(int s, int t, const vector<ivec>& enums, const vector<P>& edges, const vector<double>& ecaps, vector<int>& res)
{
	N = enums.size();
	nums.resize(enums.size());
	destNode.resize(2*edges.size());
	caps.resize(2*edges.size());
	for(int i=0; i<N; ++i) {
		nums[i].resize(enums[i].size());
		for(size_t j=0; j<enums[i].size(); ++j) {
			int e = enums[i][j];
			P p = edges[e];
			int t = i==p.first ? p.second : p.first;
			int ee = 2*e;
			if (t<i) ++ee;
			nums[i][j] = ee;
			destNode[ee] = t;
			caps[ee] = ecaps[e];
		}
	}
	heights.resize(N);
	exceeds.resize(N);

	list<int> nodes;
	heights[s] = N;
	heights[t] = 0;
	for(int i=0; i<N; ++i) {
		if (i!=s && i!=t) {
			nodes.push_back(i);
			heights[i]=1;
		}
		exceeds[i] = 0;
	}
//	cout<<"start H: "<<N<<'\n';

	for(size_t i=0; i<nums[s].size(); ++i) {
		int e = nums[s][i];
		exceeds[destNode[e]] = caps[e];
		caps[e^1] += caps[e];
		caps[e] = 0;
	}

	for(list<int>::iterator i=nodes.begin(); i!=nodes.end(); ++i) {
		int n = *i;
		if (exceeds[n]<EPS) continue;
//		cout<<"lol @ "<<n<<' '<<exceeds[n]<<' '<<heights[n]<<'\n';
		int h0 = heights[n];
		discharge(n);
		if (heights[n] != h0) {
			nodes.erase(i);
			nodes.push_front(n);
			i = nodes.begin();
		}
	}

#if 0
	cout<<"start caps:\n";
	for(size_t i=0; i<nums[s].size(); ++i) {
		int e = nums[s][i];
		cout<<destNode[e]<<' '<<caps[e]<<'\n';
	}
#endif

	used.clear();
	used.resize(N, 0);
	resDFS(s, res);
	return exceeds[t];
}

namespace {
vector<char> removed;
vector<int> iused;
vector<int> from;
vector<double> csums;
vector<vector<int> > merges;
vector<int> isEnd;
typedef pair<int,double> IDP;

void addToEdge(vector<IDP>& edges, int t, double c)
{
	size_t a=0;
	while(a<edges.size() && edges[a].first!=t) ++a;
	if (a < edges.size()) edges[a].second += c;
	else edges.push_back(IDP(t, c));
}
void addMergeTree(vector<int>& r, int n)
{
//	cout<<n<<' ';
	r.push_back(n);
	for(size_t i=0; i<merges[n].size(); ++i)
		addMergeTree(r, merges[n][i]);
}
void mergeNodes(int to, int from, vector<vector<IDP> >& edges)
{
	removed[from] = 1;
	isEnd[to] += isEnd[from];
	merges[to].push_back(from);
	for(size_t j=0; j<edges[from].size(); ++j) {
		int t = edges[from][j].first;
		if (removed[t] || t==to) continue;
		double c = edges[from][j].second;
		addToEdge(edges[to], t, c);
		addToEdge(edges[t], to, c);
	}
}
} // end anonymous namespace
double minCutFrom(int start, int start2, const vector<int>& ends, vector<vector<IDP> >& edges, vector<int>& res)
{
	typedef pair<double,int> DP;
	int N = edges.size();

	iused.clear();
	iused.resize(N, -1);
	from.clear();
	from.resize(N, -1);
	removed.clear();
	removed.resize(N, 0);
	isEnd.clear();
	isEnd.resize(N, 0);
	csums.resize(N);
	merges.resize(N);
	for(int i=0; i<N; ++i) merges[i].clear();

	for(size_t i=0; i<ends.size(); ++i) isEnd[ends[i]] = 1;

	if (start2>=0) mergeNodes(start, start2, edges);

	double bestCut=1e100;
	priority_queue<DP> pq;
	for(int i=1; i<N; ++i) {
		pq.push(DP(0, start));
		from[0]=i;
		int prev=-1, pprev=-1;
		double psum=0;
		size_t efound=0;
		while(!pq.empty()) {
			DP p = pq.top();
			pq.pop();
			int n = p.second;
			if (iused[n]==i) continue;
			iused[n] = i;
			efound += isEnd[n];
//			cout<<"juhuu "<<n<<'\n';

			for(size_t j=0; j<edges[n].size(); ++j) {
				int t = edges[n][j].first;
				if (removed[t] || iused[t]==i) continue;
				if (from[t] != i) from[t]=i, csums[t]=0;
				csums[t] += edges[n][j].second;
//				cout<<"csum "<<t<<' '<<csums[t]<<'\n';
				if (csums[t] > EPS) {
					pq.push(DP(csums[t], t));
				}
			}
			pprev = prev;
			prev = n;
			psum = p.first;
		}
		if (pprev < 0) break;
//		cout<<"lol "<<i<<' '<<prev<<' '<<pprev<<' '<<psum<<'\n';
		efound -= isEnd[prev];
		if (efound<ends.size() && psum < bestCut) {
//			cout<<"jee "<<psum<<' '<<prev<<' '<<pprev<<' '<<efound<<'/'<<ends.size()<<'\n';
			bestCut = psum;
			res.clear();
			for(int j=0; j<N; ++j)
				if (iused[j]==i && j!=prev) {
					addMergeTree(res, j);
	//				cout<<'\n';
				}

#if 0
			sort(res.begin(), res.end());
			cout<<res<<'\n';
#endif
		}
		mergeNodes(pprev, prev, edges);
	}
	return bestCut;
}
