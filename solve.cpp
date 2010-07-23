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

static int nextID=2;
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
	char* mapfile=0;
	char* dbfile=0;
	bool lprelaxation = 0;
	double maxt = 5;
	for(int i=1; i<argc; ++i) {
		char* a = argv[i];
		if (a[0]=='-') {
			switch(a[1]) {
				case 'r': lprelaxation=1; break;
				case 't': maxt = atof(argv[++i]); break;
			}
		} else if (!mapfile) mapfile = a;
		else if (!dbfile) dbfile = a;
	}
	if (mapfile) {
		ifstream in(mapfile);
		assert(in);
		readGraph(in);
	} else {
		readGraph(cin);
	}
	if (dbfile) {
		ifstream db(dbfile);
		assert(db);
		readDB(db);
	} else readDB(cin);

	vector<double> probs(purchases.size(), 1./purchases.size());

	if (lprelaxation) {
		double r = routeLP(probs);
		cout<<r<<'\n';
		return 0;
	}
	initTSPCost(purchases, probs);
	double r = antColony(expectedTotalCost, maxt);
	cout<<"Final cost "<<exactTotalCost(bestPath)<<'\n';
	cout<<bestPath<<'\n';
	cout<<r<<'\n';

//	cout<<bestPath.size()<<' ';
#if 1
	for(size_t i=0; i<bestPath.size(); ++i) {
		int k=bestPath[i];
		cout<<pos[k].x<<' '<<pos[k].y<<' ';
	}
	cout<<'\n';
#endif
}
