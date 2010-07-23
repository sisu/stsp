#include <glpk.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <cassert>
#include "tspCost.hpp"
using namespace std;

extern vector<vector<double> > dist;
extern vector<vector<int> > conn;
extern vector<vector<int> > purchases;

namespace {

typedef pair<int,int> P;
vector<double> edists;
vector<P> edges;
vector<vector<int> > enums;
int K, E, N;

glp_prob* lp;
vector<int> cols;
vector<double> row;

vector<double> vals;

const double EPS = 1e-6;

int addDegreeConstraints()
{
	int added=0;
	for(int i=0; i<=K; ++i) {
		for(int j=2; j<N; ++j) {
			double sum=0;
			for(size_t k=0; k<enums[j].size(); ++k) {
				sum += vals[enums[j][k]];
				if (i>0) sum += vals[E*i + enums[j][k]];
			}
			for(size_t k=0; k<enums[j].size(); ++k) {
				double x = vals[E*i + enums[j][k]];
				double xx = vals[enums[j][k]];
				double s = sum - (i>0 ? x+xx : x);
				if (s < x - EPS) {
//					cout<<"asd "<<i<<' '<<j<<' '<<k<<": "<<x<<' '<<sum<<' '<<s<<'\n';
					int z=0;
					for(size_t l=0; l<enums[j].size(); ++l) {
						int n = enums[j][l];
//						cout<<"ind "<<n<<'\n';
						cols[++z] = n;
						row[z] = l==k ? -1 : 1;
						if (i>0) {
							cols[++z] = E*i + n;
//							cout<<"ind' "<<E*i + n<<'\n';
							row[z] = row[z-1];
						}
					}
					int r = glp_add_rows(lp, 1);
					glp_set_row_bnds(lp, r, GLP_LO, 0, 0);
					glp_set_mat_row(lp, r, z, &cols[0], &row[0]);
					++added;
				}
			}
		}
	}
	if (added) cout<<"Added "<<added<<" degree constraints\n";
	return added;
}

int otherNode(int e, int n)
{
	P p = edges[e];
	return p.first == n ? p.second : p.first;
}

vector<int> used;
vector<double> csums;
vector<int> found;
int addSubtoursSingle(int p)
{
	used.clear();
	used.resize(N, -1);
	found.clear();
	found.resize(N, -1);
	csums.resize(N);
	vector<int> cur;

	int added=0;
	for(size_t i=0; i<purchases[p-1].size(); ++i) {
		int s = purchases[p-1][i];
		if (used[s]>=0) continue;
		cur.clear();
		found[s] = s;

		priority_queue<P> q;
		q.push(P(0,s));
		double csum=0;
		while(!q.empty()) {
			P pp = q.top();
			q.pop();
			int n = pp.second;
			if (used[n]>=0) continue;
			csum -= pp.first;
			cur.push_back(n);

			used[n] = s;
			for(size_t j=0; j<enums[n].size(); ++j) {
				int e = enums[n][j];
				int t = otherNode(e, n);
				if (used[t]>=0) {
					if (found[t]!=s)
						csum += vals[E*p + e] + vals[e];
					continue;
				}
				if (found[t]!=s) found[t]=s, csums[t] = 0;
				csum += vals[E*p + e] + vals[e];
//				cout<<"adding to csum "<<s<<' '<<csum<<' '<<vals[E*p+e]<<' '<<vals[e]<<'\n';
				if (csums[t] > .2)
					q.push(P(csums[t], t));
			}

			if (csum < 2 - EPS) {
				int z=0;
				for(size_t k=0; k<cur.size(); ++k) {
					int m = cur[k];
//					cout<<"adding vars: "<<m<<' '<<enums[m].size()<<'\n';
					for(size_t l=0; l<enums[m].size(); ++l) {
						int e = enums[m][l];
						int t = otherNode(e, m);
//						cout<<"other: "<<t<<' '<<s<<' '<<found[t]<<'\n';
						if (used[t]>=0 && found[t]==s) continue;
						cols[++z] = E*p + e;
						row[z] = 1;
						cols[++z] = e;
						row[z] = 1;
					}
				}
//				cout<<"adding subtour cut "<<p<<' '<<cur.size()<<' '<<s<<' '<<z<<' '<<csum<<'\n';
				double tsum=0;
				for(int i=1; i<=z; ++i) tsum += vals[cols[i]];
				cout<<tsum<<'\n';

				assert(z);
				int r = glp_add_rows(lp, 1);
				glp_set_row_bnds(lp, r, GLP_LO, 2, 0);
				glp_set_mat_row(lp, r, z, &cols[0], &row[0]);
				++added;
			}
		}
	}
	if (added) cout<<"added "<<added<<" subtour cuts from "<<p<<'\n';
	return added;
}
int addSubtourConstraints()
{
	int added = 0;
	for(int i=1; i<=K; ++i) {
		added += addSubtoursSingle(i);
	}
	return added;
}

bool addConstraints()
{
	if (addDegreeConstraints()) return 1;
	if (addSubtourConstraints()) return 1;
	return 0;
}

} // end anonymous namespace

double routeLP(const vector<double>& probs)
{
	K = purchases.size();
	E = 0;
	N = conn.size();;
	enums.resize(1+N);
	edists.resize(1); // dummy to start real indexing from 1
	edges.resize(1);
	for(size_t i=0; i<conn.size(); ++i) {
		for(size_t j=0; j<conn[i].size(); ++j) {
			size_t t = conn[i][j];
			if (t < i) continue;
			int n = edges.size();
			enums[i].push_back(n);
			enums[t].push_back(n);
			edges.push_back(P(i,t));
//			cout<<"setting edists "<<edists.size()<<": "<<dist[i][t]<<'\n';
			edists.push_back(dist[i][t]);
		}
	}
	E = edges.size() - 1;
	lp = glp_create_prob();
	int vs = E*(K+1);
	glp_add_cols(lp, vs);
	for(int i=1; i<=vs; ++i) {
		glp_set_col_bnds(lp, i, GLP_DB, 0, 1);
	}
	cols.resize(vs+1);
	row.clear();
	row.resize(vs+1, 1);

	glp_set_obj_dir(lp, GLP_MIN);
	for(int i=1; i<=E; ++i)
		glp_set_obj_coef(lp, i, edists[i] * LENGTH_FACTOR);
	for(int i=1; i<=K; ++i) {
		for(int j=1; j<=E; ++j) {
			glp_set_obj_coef(lp, E*i + j, edists[j] * probs[i-1]);
//			cout<<"setting coeff "<<i<<' '<<j<<": "<<edists[j]*probs[i-1]<<' '<<edists[j]<<'\n';
		}
	}

	{
		// Bounds for start and end nodes
		int r = glp_add_rows(lp, 2);
		for(int i=0; i<2; ++i, ++r) {
			for(size_t k=0; k<enums[i].size(); ++k) {
				cols[1+k] = enums[i][k];
			}
			glp_set_row_bnds(lp, r, GLP_FX, 1, 1);
			glp_set_mat_row(lp, r, enums[i].size(), &cols[0], &row[0]);
		}
	}

	cout<<"vars: "<<vs<<" ; "<<E<<' '<<K<<'\n';

	for(int i=0; i<K; ++i) {
		int r = glp_add_rows(lp, purchases[i].size());
		for(size_t j=0; j<purchases[i].size(); ++j, ++r) {
			int t = purchases[i][j];
//			cout<<"lol "<<i<<' '<<j<<' '<<t<<' '<<enums[t].size()<<'\n';
			int z=0;
			for(size_t k=0; k<enums[t].size(); ++k) {
				cols[++z] = enums[t][k];
				cols[++z] = enums[t][k] + E*(1+i);
			}
			glp_set_row_bnds(lp, r, GLP_LO, 2, 0);
			glp_set_mat_row(lp, r, z, &cols[0], &row[0]);
		}
	}

	glp_smcp parm;
	glp_init_smcp(&parm);
//	parm.msg_lev = GLP_MSG_ERR;
	parm.meth = GLP_PRIMAL;
	parm.presolve = 1;

	vals.resize(1 + vs);
	do {
		int r = glp_simplex(lp, &parm);
		if (r || glp_get_status(lp)!=GLP_OPT) {
			cout<<"FAILED SOLVING ROUTE RELAXATION: "<<r<<' '<<glp_get_status(lp)<<'\n';
			abort();
		}
		parm.meth = GLP_DUAL;
		parm.presolve = 0;

		cout<<"tmp result "<<glp_get_obj_val(lp)<<'\n';

		for(int i=1; i<=vs; ++i)
			vals[i] = glp_get_col_prim(lp, i);
	} while(addConstraints());

	double res = glp_get_obj_val(lp);

	glp_delete_prob(lp);

	return res;
}
