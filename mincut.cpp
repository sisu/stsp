#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
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

};

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
