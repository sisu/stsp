#include <vector>
#include <iostream>
#include <tr1/unordered_map>
#include <queue>
#include "util.hpp"
#include "Vector.hpp"
#include "Rect.hpp"
using namespace std;

extern vector<vector<int> > conn;
extern vector<vector<int> > purchases;
extern vector<Vec2> pos;
extern vector<int> itemID;

typedef pair<int,int> IP;
typedef tr1::unordered_map<IP,double> EdgeDist;
vector<EdgeDist> edgeDist;

vector<double> itemFreq;

double pathLength(const vector<int>& path)
{
	double r=0;
	for(size_t i=1; i<path.size(); ++i)
		r += length(pos[path[i]]-pos[path[i-1]]);
	return r;
}
double pathDist(int c, const vector<int>& path)
{
	double d=1e100;
	for(size_t i=1; i<path.size(); ++i)
		d = min(d, edgeDist[c][IP(path[i-1],path[i])]);
	return d;
}

double distanceCost(const vector<int>& path)
{
//	cout<<path<<'\n';
//	return path.size();

	const double TOTAL_DIST_FACT = .5;

	double r = TOTAL_DIST_FACT * pow(pathLength(path),1);

	for(size_t i=0; i<itemFreq.size(); ++i) {
		r += pathDist(i,path)*itemFreq[i];
//		r += pow(pathDist(i,path),1+itemFreq[i]);
	}
	return r;
}
double distanceCost2(const vector<int>& path)
{
	typedef pair<double,int> P;
	vector<P> v;
	for(size_t i=0; i<itemFreq.size(); ++i) {
		double p = itemFreq[i] * pathDist(i,path);
		v.push_back(P(p,i));
	}
	sort(v.begin(),v.end());
	double r = .5 * pathLength(path);
	double a=2;
	for(size_t i=0; i<v.size(); ++i, a*=.98) {
		r += v[i].first * a;
	}
	return r;
}

vector<vector<pair<IP,double> > > straightEdges;

EdgeDist getEdgeDists(int s)
{
	EdgeDist res;
	typedef pair<double,int> P;
	priority_queue<P,vector<P>,greater<P> > q;
	q.push(P(0,s));

	vector<char> used(pos.size());
	while(!q.empty()) {
		P p=q.top();
		q.pop();
		int n=p.second;
		if (used[n]) continue;
		used[n]=1;

		double d = p.first;
		for(size_t i=0; i<conn[n].size(); ++i) {
			int t = conn[n][i];
			if (used[t]) continue;
			IP pp(n,t);
			if (!res.count(pp) || res[pp]>d) res[pp]=res[IP(t,n)]=d;
			q.push(P(d+length(pos[t]-pos[n]), t));
		}
#if 1
		for(size_t i=0; i<straightEdges[n].size(); ++i) {
			IP e = straightEdges[n][i].first;
			double dd = d + straightEdges[n][i].second;
			if (!res.count(e) || res[e]>dd) res[e]=res[IP(e.second,e.first)]=dd;
		}
#endif
	}
	return res;
}
extern vector<Rect> bigRects;
void genStraightEdges()
{
	// FIXME: O(n^4)...
	straightEdges.resize(pos.size());
	for(size_t i=0; i<pos.size(); ++i) {
		for(size_t j=0; j<pos.size(); ++j) {
			if (i==j) continue;
			for(size_t k=0; k<conn[j].size(); ++k) {
				size_t t=conn[j][k];
				if (t<j || t==i) continue;
				if (!between(pos[j],pos[t],pos[i])) continue;
				Vec2 p = projectionToLine(pos[j],pos[t],pos[i]);
				bool ok=1;
				for(size_t l=0; l<bigRects.size(); ++l) {
					Rect r = bigRects[l];
					if (hitsRect(pos[i],p,r)) {
						ok=0;
						break;
					}
				}
				if (ok) {
//					cout<<"straight edge "<<i<<' '<<j<<' '<<t<<'\n';
					double d = length(p-pos[i]);
					straightEdges[i].push_back(make_pair(IP(j,t), d));
				}
			}
		}
	}
}

void initDistances()
{
	genStraightEdges();
	edgeDist.resize(itemID.size());
	for(size_t i=0; i<itemID.size(); ++i) {
		edgeDist[i] = getEdgeDists(itemID[i]);
	}

	itemFreq.resize(itemID.size());
	for(size_t i=0; i<purchases.size(); ++i) {
		for(size_t j=0; j<purchases[i].size(); ++j)
			++itemFreq[purchases[i][j]];
	}
	for(size_t i=0; i<itemFreq.size(); ++i)
		itemFreq[i] /= purchases.size();
}
