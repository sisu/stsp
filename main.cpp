#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <cctype>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "Vector.hpp"
#include "Rect.hpp"
#include "distanceCost.hpp"
#include "util.hpp"
using namespace std;

double antColony(double(*)(const vector<int>&));


double W,H;

vector<vector<int> > conn;
vector<Vec2> pos;
int startI, endI;

vector<char> used;
vector<char> onPath;

vector<vector<int> > purchases;

vector<int> bestPath;

struct Commodity {
	int id;
	Vec2 pos;
};
vector<Commodity> goods;
vector<int> itemID;

static int nextID=0;
static map<string,int> IDMap;
int getID(const string s)
{
	if (IDMap.count(s)) return IDMap[s];
	return IDMap[s]=nextID++;
}
int tryGetID(const string s)
{
	if (IDMap.count(s)) return IDMap[s];
	return -1;
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

vector<Rect> rects;

const double CX=4, CY=4;
bool outside(Vec2 v)
{
	return v.x<CX || v.y<CY || v.x>W-CX || v.y>H-CY;
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
Vec2 startV,endV;
vector<Rect> bigRects;
void genGraph()
{
	bigRects = rects;
	for(size_t i=0; i<rects.size(); ++i) {
		Rect& b = bigRects[i];
		Rect r = rects[i];
		b.x1 = r.x1-CX;
		b.x2 = r.x2+CX;
		b.y1 = r.y1-CY;
		b.y2 = r.y2+CY;
	}

//	vector<Segment> walls;
	conn.resize(4*rects.size() + goods.size() + 2);
	pos.resize(conn.size());
	for(size_t i=0; i<rects.size(); ++i) {
		genCorners(bigRects[i], (&pos[0])+4*i);
	}
	size_t s0 = 4*rects.size();
	itemID.resize(nextID);
	for(size_t i=0; i<goods.size(); ++i) {
		Commodity c = goods[i];
		pos[s0+i] = c.pos;
		itemID[i] = s0+i;
	}
	startI = s0+goods.size();
	endI = startI+1;
	pos[startI] = startV;
	pos[endI] = endV;

	for(size_t i=0; i<pos.size(); ++i)
		for(size_t j=0; j<i; ++j)
			checkedJoin(i, j, bigRects);

	for(size_t i=0; i<conn.size(); ++i)
		sort(conn[i].begin(),conn[i].end());
}

void readArea(istream& in)
{
	in>>W>>H;
	in>>startV.x>>startV.y>>endV.x>>endV.y;

#if 0
	string s;
	while(in>>s) {
		if (isdigit(s[0])) {
			Rect r;
			r.x1 = atof(s.c_str());
			in>>r.y1>>r.x2>>r.y2;
			rects.push_back(r);
		} else {
			Commodity c;
			c.id = getID(s);
			in>>c.pos.x>>c.pos.y;
			goods.push_back(c);
		}
	}
#endif
	int n;
	in>>n;
	rects.resize(n);
	for(int i=0; i<n; ++i) {
		Rect& r=rects[i];
		in>>r.x1>>r.y1>>r.x2>>r.y2;
	}
	int k;
	in>>k;
	goods.resize(k);
	string s;
	for(int i=0; i<k; ++i) {
		in>>s;
		Commodity& c = goods[i];
		c.id = getID(s);
		in>>c.pos.x>>c.pos.y;
	}
}
void readDB(istream& in)
{
#if 0
	string s;
	while(getline(in,s)) {
		if (s.empty()) continue;
		istringstream ss(s);
		vector<int> v;
		string a;
		while(ss>>a) v.push_back(getID(a));
		purchases.push_back(v);
	}
#endif
	int n;
	in>>n;
	purchases.resize(n);
	for(int i=0; i<n; ++i) {
		int k;
		in>>k;
		purchases[i].reserve(k);
		string s;
		for(int j=0; j<k; ++j) {
			in>>s;
			int id = tryGetID(s);
			if (id>=0)
				purchases[i].push_back(id);
		}
	}

}

int main(int argc, char* argv[])
{
	srand(time(0));
	if (argc>2) {
		ifstream in(argv[1]);
		readArea(in);
		ifstream db(argv[2]);
		readDB(db);
	} else {
		readArea(cin);
		readDB(cin);
	}
	genGraph();

	initDistances();

//	double r = optimize(distanceCost);
//	double r = antColony(distanceCost);
//	double r = antColony(distanceCost2);
	double r = antColony(tspCost);
	cout<<r<<'\n';
	cout<<bestPath<<'\n';

//	cout<<bestPath.size()<<' ';
	for(size_t i=0; i<bestPath.size(); ++i) {
		int k=bestPath[i];
		cout<<pos[k].x<<' '<<pos[k].y<<' ';
	}
	cout<<'\n';
}
