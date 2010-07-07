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
#include "Area.hpp"
using namespace std;

double antColony(double(*)(const vector<int>&), double t);


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

//	vector<Segment> walls;
	conn.resize(4*rects.size() + items.size() + 2);
	pos.resize(conn.size());
	for(size_t i=0; i<rects.size(); ++i) {
		genCorners(bigRects[i], (&pos[0])+4*i);
	}
	size_t s0 = 4*rects.size();
	itemID.resize(items.size());
	int cur=0;
	for(map<string,Vec2>::iterator i=items.begin(); i!=items.end(); ++i, ++cur) {
		int x = s0+cur;
		pos[x] = i->second;
		itemID[cur] = x;
	}
	/*
	for(size_t i=0; i<goods.size(); ++i) {
		Commodity c = goods[i];
		pos[s0+i] = c.pos;
		itemID[i] = s0+i;
	}*/

	startI = s0+items.size();
	endI = startI+1;
	pos[startI] = area.startV;
	pos[endI] = area.endV;

	for(size_t i=0; i<pos.size(); ++i)
		for(size_t j=0; j<i; ++j)
			checkedJoin(i, j, bigRects);

	for(size_t i=0; i<conn.size(); ++i)
		sort(conn[i].begin(),conn[i].end());
}

void readDB(istream& in)
{
	string s;
	while(getline(in,s)) {
		if (s.empty()) continue;
		istringstream ss(s);
		vector<int> v;
		string a;
		while(ss>>a) {
			int x = tryGetID(a);
			if (x>=0)
				v.push_back(x);
		}
		purchases.push_back(v);
	}
}
void readInput(istream& in)
{
	area.readSimple(in);
	for(map<string,Vec2>::iterator i=area.items.begin(); i!=area.items.end(); ++i)
		getID(i->first);

	int m;
	in>>m;
	purchases.resize(m);
	for(int i=0; i<m; ++i) {
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
		assert(in);
//		readArea(in);
		area.read(in);
		for(map<string,Vec2>::iterator i=area.items.begin(); i!=area.items.end(); ++i)
			getID(i->first);

		ifstream db(argv[2]);
		assert(db);
		readDB(db);
	} else {
		readInput(cin);
	}
	genGraph();

	cout<<"Graph generated "<<conn.size()<<'\n';;

	initDistances();

	cout<<"Preprocessing done\n";

//	double r = antColony(distanceCost);
//	double r = antColony(distanceCost2);
	double r = antColony(distanceCost3, 3);
//	double r = antColony(tspCost);
	cout<<r<<'\n';
	cout<<bestPath<<'\n';

//	cout<<bestPath.size()<<' ';
	for(size_t i=0; i<bestPath.size(); ++i) {
		int k=bestPath[i];
		cout<<pos[k].x<<' '<<pos[k].y<<' ';
	}
	cout<<'\n';
}
