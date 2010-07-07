#include <vector>
#include <iostream>
#include <tr1/unordered_map>
#include <queue>
#include <algorithm>
#include "util.hpp"
#include "Vector.hpp"
#include "Rect.hpp"
#include "tsp.hpp"
using namespace std;

typedef vector<int> ivec;
vector<ivec> conn;
vector<vector<double> > edist;
ivec itemID;
const int startI=0, endI=1;

typedef pair<int,int> IP;
//typedef tr1::unordered_map<IP,double> EdgeDist;
typedef vector<double> EdgeDist;
tr1::unordered_map<IP,int> edgeNums;
vector<EdgeDist> edgeDist;

void genEdgeNums()
{
	int n=0;
	for(size_t i=0; i<conn.size(); ++i) {
		for(size_t j=0; j<conn[i].size(); ++j) {
			int t=conn[i][j];
			if (edgeNums.count(IP(i,t))) continue;
			edgeNums[IP(i,t)] = edgeNums[IP(t,i)] = n++;
		}
	}
	cout<<"Edges: "<<n<<'\n';
}
vector<vector<double> > dist;
vector<ivec> to;

int main()
{
	int N,K;
	cin>>N>>K;
	cout<<"asd "<<itemID<<'\n';
	genEdgeNums();

	for(int i=0; i<N; ++i)
		for(int j=0; j<N; ++j)
			for(int k=0; k<N; ++k)
				if (dist[j][i] + dist[i][k] < dist[j][k]) {
					dist[j][k] = dist[j][i] + dist[i][k];
					to[j][k] = to[j][i];
				}
}
