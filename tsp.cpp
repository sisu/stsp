#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <queue>
#include <iostream>
#include "Vector.hpp"
#include "util.hpp"
using namespace std;

extern vector<Vec2> pos;
extern vector<vector<int> > conn;
extern vector<int> itemID;
extern int startI, endI;

typedef vector<int> ivec;

vector<vector<vector<int> > > itemPath;
vector<vector<int> > startPath, endPath;

namespace {

vector<vector<double> > itemDist;
vector<double> startDist, endDist;
double pathCost(const ivec& v)
{
	double r = startDist[v[0]] + endDist[v.back()];
	for(size_t i=1; i<v.size(); ++i)
		r += itemDist[v[i-1]][v[i]];
	return r;
}

double randf()
{
	return rand()/(double)RAND_MAX;
}

vector<int> rItemID;

vector<bool> used;
vector<int> from;

void genShortestPaths(int sn, vector<double>& dists, vector<ivec>& paths)
{
	int s = itemID[sn];

	used.resize(0);
	used.resize(conn.size(),0);
	from.resize(conn.size());

	typedef pair<int,int> IP;
	typedef pair<double,IP> P;
	priority_queue<P,vector<P>,greater<P> > q;
	q.push(P(0,IP(s,-1)));
	while(!q.empty()) {
		P p = q.top();
		q.pop();
		double d=p.first;
		int n=p.second.first;
		int f=p.second.second;
		if (used[n]) continue;
		used[n]=1;
		from[n]=f;
		if (rItemID[n]>=0 || n==startI || n==endI) {
			ivec& path = n==startI ? startPath[sn] : n==endI ? endPath[sn] : paths[rItemID[n]];
			int nn = n;
			do {
//				cout<<"lol "<<sn<<' '<<n<<' '<<rItemID[n]<<' '<<nn<<'\n';
				path.push_back(nn);
			} while((nn = from[nn])>=0);
			if (n!=startI) reverse(path.begin(),path.end());
//			if (n==startI) cout<<"found start path for "<<sn<<" : "<<path<<'\n';

			if (n==startI) startDist[sn]=d;
			else if (n==endI) endDist[sn]=d;
			else dists[rItemID[n]] = d;
		}

		for(size_t i=0; i<conn[n].size(); ++i) {
			int x=conn[n][i];
			if (used[x]) continue;
			q.push(P(d+length(pos[x]-pos[n]), IP(x,n)));
		}
	}
}

}

ivec TSP(const ivec& purchases)
{
	if (purchases.size()<2) return ivec(purchases);
	ivec cur = purchases;
	ivec best = cur;
	double bdist = pathCost(cur);
	double cdist = bdist;
	for(double t=10; t>.1; t*=.999) {
		size_t a = rand()%(cur.size()-1);
		size_t b = rand()%(cur.size()-1);
		if (b>=a) ++b;
		else swap(a,b);

		int ca=cur[a], cb=cur[b];
		double da = a ? itemDist[cur[a-1]][cb] -  itemDist[cur[a-1]][ca] : startDist[cb] - startDist[ca];
		double db = b<cur.size()-1 ? itemDist[ca][cur[b+1]] - itemDist[cb][cur[b+1]] : endDist[ca] - endDist[cb];
		double d = da+db;
		if (d<0 || randf()<exp(-d/t)) {
			reverse(cur.begin()+a,cur.begin()+b+1);
			cdist = pathCost(cur);
			if (cdist < bdist) {
				best = cur;
				bdist = cdist;
			}
		}
	}
	return best;
}

void initTSP()
{
	rItemID.resize(pos.size(), -1);
	for(size_t i=0; i<itemID.size(); ++i)
		rItemID[itemID[i]] = i;

	itemDist.resize(itemID.size());
	itemPath.resize(itemID.size());
	startDist.resize(itemID.size());
	startPath.resize(itemID.size());
	endDist.resize(itemID.size());
	endPath.resize(itemID.size());
	for(size_t i=0; i<itemDist.size(); ++i) {
		itemDist[i].resize(itemDist.size());
		itemPath[i].resize(itemDist.size());
		genShortestPaths(i, itemDist[i], itemPath[i]);
	}
}
