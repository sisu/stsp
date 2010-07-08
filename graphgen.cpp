#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <cmath>
#include "Vector.hpp"
#include "Rect.hpp"
#include "Area.hpp"
using namespace std;

Area area;

vector<vector<int> > conn;
vector<Vec2> pos;
int startI, endI;

vector<char> used;
vector<char> onPath;

vector<vector<int> > purchases;

vector<int> bestPath;

vector<int> itemID;

static int nextID=0;
static map<string,int> IDMap;
int getID(const string s)
{
	if (IDMap.count(s)) return IDMap[s];
	return IDMap[s]=nextID++;
}

void genCorners(Rect r, Vec2 vs[4])
{
	double x1=r.x1, x2=r.x2;
	double y1=r.y1, y2=r.y2;
	vs[0]=Vec2(x1,y1);
	vs[1]=Vec2(x2,y1);
	vs[2]=Vec2(x2,y2);
	vs[3]=Vec2(x1,y2);
}

const double CX=4, CY=4;
bool outside(Vec2 v)
{
	return v.x<CX || v.y<CY || v.x>area.W-CX || v.y>area.H-CY;
}
bool hitsAnyRect(Vec2 a, Vec2 b, const vector<Rect>& rects)
{
	if (outside(a) || outside(b)) return 1;
	for(size_t i=0; i<rects.size(); ++i) {
		if (hitsRect(a,b,rects[i])) return 1;
	}
	return 0;
}
void checkedJoin(int a, int b, const vector<Rect>& rects)
{
	if (!hitsAnyRect(pos[a], pos[b], rects)) {
		conn[a].push_back(b);
		conn[b].push_back(a);
//		cout<<"conn "<<a<<' '<<b<<" ; "<<pos[a]<<' '<<pos[b]<<'\n';
	}
}
vector<Rect> bigRects;
void genGraph()
{
	vector<Rect>& rects = area.rects;
	map<string,Vec2>& items = area.items;
	bigRects = rects;
	for(size_t i=0; i<rects.size(); ++i) {
		Rect& b = bigRects[i];
		Rect r = rects[i];
		b.x1 = r.x1-CX;
		b.x2 = r.x2+CX;
		b.y1 = r.y1-CY;
		b.y2 = r.y2+CY;
	}

	pos.resize(2+items.size());
	itemID.resize(items.size());
	pos[startI = 0] = area.startV;
	pos[endI = 1] = area.endV;
	int cur=0;
	for(map<string,Vec2>::iterator i=items.begin(); i!=items.end(); ++i, ++cur) {
		int x = 2+cur;
		pos[x] = i->second;
		itemID[cur] = x;
	}
	size_t s0 = pos.size();

//	vector<Segment> walls;
	pos.resize(s0 + 4*rects.size());
	for(size_t i=0; i<rects.size(); ++i) {
		genCorners(bigRects[i], &pos[s0 + 4*i]);
	}
//	cout<<"before "<<pos.size()<<'\n';
	pos.erase(remove_if(pos.begin()+s0, pos.end(), outside), pos.end());
//	cout<<"after "<<pos.size()<<'\n';

	conn.resize(pos.size());

	for(size_t i=0; i<pos.size(); ++i)
		for(size_t j=0; j<i; ++j)
			checkedJoin(i, j, bigRects);

	for(size_t i=0; i<conn.size(); ++i)
		sort(conn[i].begin(),conn[i].end());
}

void genDists()
{
	int N = conn.size();
	vector<vector<double> > dist(N, vector<double>(N, 1e50));
	vector<vector<int> > to(N, vector<int>(N));
	for(int i=0; i<N; ++i) {
		dist[i][i] = 0;
		to[i][i] = i;
		for(size_t j=0; j<conn[i].size(); ++j) {
			int t = conn[i][j];
			dist[i][t] = length(pos[i] - pos[t]);
			to[i][t] = t;
		}
	}
	for(int i=0; i<N; ++i)
		for(int j=0; j<N; ++j)
			for(int k=0; k<N; ++k)
				if (dist[j][i] + dist[i][k] < dist[j][k]) {
					dist[j][k] = dist[j][i] + dist[i][k];
					to[j][k] = to[j][i];
				}
	for(int i=0; i<N; ++i) {
		for(int j=0; j<N; ++j)
			cout<<dist[i][j]<<' '<<to[i][j]<<' ';
		cout<<'\n';
	}
}

int main(int argc, char* argv[])
{
	area.read(cin);
	for(map<string,Vec2>::iterator i=area.items.begin(); i!=area.items.end(); ++i)
		getID(i->first);
	genGraph();

	cout<<conn.size()<<' '<<area.items.size()<<'\n';
	for(map<string,Vec2>::iterator i=area.items.begin(); i!=area.items.end(); ++i)
		cout<<i->first<<'\n';
	for(size_t i=0; i<conn.size(); ++i) {
		cout<<conn[i].size()<<' '<<pos[i].x<<' '<<pos[i].y<<'\n';
		for(size_t j=0; j<conn[i].size(); ++j)
			cout<<conn[i][j]<<' '<<length(pos[i]-pos[conn[i][j]])<<' ';
		cout<<'\n';
	}
	genDists();
}
