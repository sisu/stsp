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

double W,H;

vector<vector<int> > conn;
vector<Vec2> pos;
int startI, endI;

vector<char> used;
vector<char> onPath;

vector<vector<int> > purchases;

double antColony(double(*)(const vector<int>&));

bool startDFS(int s, int t, vector<int>& out)
{
	used[s]=1;
	if (s==t) {
		out.push_back(s);
		return 1;
	}
	for(size_t i=0; i<conn[s].size(); ++i) {
		int x = conn[s][i];
		if (used[x]) continue;
		if (startDFS(x, t, out)) {
			out.push_back(s);
			return 1;
		}
	}
	return 0;
}
bool augmentDFS(int s, vector<int>& out, int no1, int no2)
{
	used[s]=1;
	if (no1<0 && no2<0 && onPath[s]) {
		out.push_back(s);
		return 1;
	}
	random_shuffle(conn[s].begin(),conn[s].end());
	for(size_t i=0; i<conn[s].size(); ++i) {
//		swap(conn[s][i],conn[s][i+rand()%(conn[s].size()-i)]);
		int x = conn[s][i];
		if (x==no1 || x==no2) continue;
		if (used[x]) continue;
		if (augmentDFS(x, out, -1, -1)) {
			out.push_back(s);
			return 1;
		}
	}
	return 0;
}

double randf()
{
	return rand()/(double)RAND_MAX;
}

vector<int> bestPath;
template<class F>
double optimize(F cost)
{
	used.resize(pos.size());
	vector<int> path;
	assert( startDFS(startI, endI, path) );
	reverse(path.begin(),path.end());

	bestPath = path;
	double best = cost(path);

	onPath.resize(pos.size());
	for(size_t i=0; i<path.size(); ++i)
		onPath[path[i]] = 1;

	vector<int> augment;
	vector<int> tmpV;
	double curC = best;

	for(double t=10; t>.1; t*=.9999) {
//		cout<<"cur "<<path<<'\n';

		fill(used.begin(),used.end(),0);
		augment.clear();
		int s = rand() % (path.size()-1);
		if (!augmentDFS(path[s], augment, path[s+1], s>0?path[s-1]:-1)) continue;

		int e=0;
		while(path[e]!=augment[0]) ++e;
//		cout<<"augment "<<augment<<" ; "<<s<<' '<<e<<'\n';
		if (e<s)
			swap(s,e);
		else
			reverse(augment.begin(),augment.end());

		tmpV.clear();
		tmpV.insert(tmpV.end(), path.begin(), path.begin()+s);
		tmpV.insert(tmpV.end(), augment.begin(), augment.end());
		tmpV.insert(tmpV.end(), path.begin()+e+1, path.end());

		double c = cost(tmpV);
		if (c < curC || randf() < exp((c-curC)/t)) {
			path.swap(tmpV);
			curC = c;

			if (curC < best) {
				best=curC;
				bestPath = path;
				cout<<"new best "<<best<<" : "<<bestPath<<'\n';
//				cout<<s<<' '<<e<<'\n';

#if 0
				for(size_t i=1; i<path.size(); ++i) {
					int a=path[i-1], b=path[i];
//					cout<<"trying "<<a<<' '<<b<<'\n';
					assert(find(conn[a].begin(),conn[a].end(),b)!=conn[a].end());
				}
#endif
			}

			fill(onPath.begin(),onPath.end(),0);
			for(size_t i=0; i<path.size(); ++i)
				onPath[path[i]] = 1;
		}
	}
	return best;
}

struct Commodity {
	int id;
	Vec2 pos;
};
vector<Commodity> goods;
vector<int> itemID;

static int nextID=0;
int getID(const string s)
{
	static map<string,int> ID;
	if (ID.count(s)) return ID[s];
	return ID[s]=nextID++;
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
		purchases[i].resize(k);
		string s;
		for(int j=0; j<k; ++j) {
			in>>s;
			purchases[i][j] = getID(s);
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
	double r = antColony(distanceCost);
//	double r = antColony(distanceCost2);
	cout<<r<<'\n';
	cout<<bestPath<<'\n';

//	cout<<bestPath.size()<<' ';
	for(size_t i=0; i<bestPath.size(); ++i) {
		int k=bestPath[i];
		cout<<pos[k].x<<' '<<pos[k].y<<' ';
	}
	cout<<'\n';
}
