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
#include "util.hpp"
#include "tspCost.hpp"
#include "glproute.hpp"
using namespace std;

double antColony(double(*)(const vector<int>&), double time);

vector<vector<int> > conn;
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

void readDB(istream& in)
{
	string s;
	while(getline(in,s)) {
		if (s.empty()) continue;
		istringstream ss(s);
		vector<int> v;
		string a;
		while(ss>>a) {
			if (a==";") return;
			int x = tryGetID(a);
			if (x>=0)
				v.push_back(x);
		}
		if (v.empty()) continue;
		purchases.push_back(v);
	}
}
vector<vector<int> > nextP;
vector<vector<double> > dist;
vector<Vec2> pos;
void readGraph(istream& in)
{
	int N,K;
	in>>N>>K;
	string name;
	for(int i=0; i<K; ++i) {
		in>>name;
		getID(name);
	}
	conn.resize(N);
	pos.resize(N);
	for(int i=0; i<N; ++i) {
		int k;
		in>>k>>pos[i].x>>pos[i].y;
		for(int j=0; j<k; ++j) {
			int t;
			double d;
			in>>t>>d;
			conn[i].push_back(t);
		}
	}
	nextP.resize(N, vector<int>(N));
	dist.resize(N, vector<double>(N));
	for(int i=0; i<N; ++i)
		for(int j=0; j<N; ++j)
			in>>dist[i][j]>>nextP[i][j];
	startI = 0;
	endI = 1;
}
#if 0
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
#endif

int main(int argc, char* argv[])
{
	srand(time(0));
	if (argc>1) {
		ifstream in(argv[1]);
		assert(in);
		readGraph(in);
	} else {
		readGraph(cin);
	}
	if (argc > 2) {
		ifstream db(argv[2]);
		assert(db);
		readDB(db);
	} else readDB(cin);

	vector<double> probs(purchases.size(), 1./purchases.size());
#if 1
	initTSPCost(purchases, probs);
	double r = antColony(expectedTotalCost, 5);
#else
	double r = routeLP(probs);
#endif
	cout<<"Final cost "<<exactTotalCost(bestPath)<<'\n';
	cout<<r<<'\n';
	cout<<bestPath<<'\n';

//	cout<<bestPath.size()<<' ';
#if 1
	for(size_t i=0; i<bestPath.size(); ++i) {
		int k=bestPath[i];
		cout<<pos[k].x<<' '<<pos[k].y<<' ';
	}
	cout<<'\n';
#endif
}
