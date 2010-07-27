#include <glpk.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <cassert>
#include <cmath>
#include <iomanip>
#include "tspCost.hpp"
#include "util.hpp"
using namespace std;

extern vector<vector<double> > edgeDist;
extern vector<vector<int> > conn;
extern vector<vector<int> > purchases;
extern bool robustOpt;
extern bool singleDir;
extern bool singleDirAll;

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
	if (!singleDir) return 0;
	int added=0;
	int end = singleDirAll ? K : 0;
	for(int i=0; i<=end; ++i) {
		for(int j=2; j<N; ++j) {
			double sum=0;
			for(size_t k=0; k<enums[j].size(); ++k) {
				int e = enums[j][k];
				sum += vals[e] + (i>0 ? vals[E*i + e] : 0);
			}
			if (sum > 2+EPS) {
				int z=0;
				for(size_t l=0; l<enums[j].size(); ++l) {
					int n = enums[j][l];
					cols[++z] = n;
					row[z] = 1;
					if (i>0) {
						cols[++z] = E*i+n;
						row[z] = 1;
					}
				}
				int r = glp_add_rows(lp, 1);
				glp_set_row_bnds(lp, r, GLP_UP, 0, 2);
				glp_set_mat_row(lp, r, z, &cols[0], &row[0]);
				++added;
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
//					cout<<"adding deg constr "<<i<<' '<<j<<' '<<z<<' '<<k<<' '<<sum<<'\n';
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

void dumpGraph(int p=0)
{
	int pp = cout.precision();
	for(int i=1; i<=E; ++i) {
		double v = vals[p*E + i];
		if (p) v+=vals[i];
		if (v<EPS) continue;
		P p = edges[i];
		cout<<p.first<<' '<<p.second<<" "<<setprecision(3)<<v<<'\n';
	}
//	cout<<resetiosflags(ios_base::precision);
	cout.precision(pp);
}

vector<int> used;
vector<double> csums;
vector<int> found;
int addSubtoursSingle(int p, bool weakCuts=0)
{
	used.clear();
	used.resize(N, -1);
	found.clear();
	found.resize(N, -1);
	csums.resize(N);
	vector<int> cur;

	typedef pair<double,int> DP;

	int added=0;
	for(size_t i=0; (p>0 && i<purchases[p-1].size()) || (p==0 && (int)i<N); ++i) {
		int s = p>0 ? purchases[p-1][i] : i;
		if (used[s]>=0) continue;
		if (p==0 && i>1) {
			double csum=0;
			for(size_t j=0; j<enums[i].size(); ++j) csum += vals[enums[i][j]];
			if (csum < EPS) continue;
		}
		cur.clear();
		found[s] = s;

//		if (p==0 && s>1) cout<<"Trying to find main path subtour "<<s<<'\n';

		int rownum = -1;
		double maxViol = -1;

		priority_queue<DP> q;
		q.push(DP(0,s));
		double csum=0;
		double isum=0;
		while(!q.empty()) {
			DP pp = q.top();
			q.pop();
			int n = pp.second;
			if (used[n]>=0) continue;
//			cout<<"popping "<<n<<' '<<pp.first<<' '<<csum<<' '<<isum<<'\n';
			csum -= pp.first;
			isum += pp.first;
			cur.push_back(n);

			used[n] = s;
			if (n!=s && (n==0 || n==1)) break;
//			if (p==0 && s==0) cout<<"lol "<<s<<' '<<n<<' '<<csum<<' '<<isum<<' '<<pp.first<<'\n';
			for(size_t j=0; j<enums[n].size(); ++j) {
				int e = enums[n][j];
				int t = otherNode(e, n);
				double v = p==0 ? vals[e] : vals[e] + vals[E*p+e];
				if (used[t]>=0) {
					if (found[t]!=s)
						csum += v;
					continue;
				}
				if (found[t]!=s) found[t]=s, csums[t] = v;
				else csums[t] += v;
				csum += v;
//				cout<<"adding to csum "<<s<<' '<<csum<<' '<<vals[E*p+e]<<' '<<vals[e]<<'\n';
				if (csums[t] > EPS) {
//				if (csums[t] > 0) {
					q.push(DP(csums[t], t));
//					cout<<"pushing "<<t<<' '<<csums[t]<<'\n';
				}
			}

			if ((p>0 && csum < 2 - EPS) || (p==0 && i<2 && csum < 1-EPS)) {
				int z=0;
				if (p==0) cout<<"adding start/end constr "<<csum<<' '<<isum<<'\n';
//				for(size_t k=0; k<cur.size(); ++k) used[cur[k]]=found[cur[k]]=s;
				for(size_t k=0; k<cur.size(); ++k) {
					int m = cur[k];
//					cout<<"adding vars: "<<m<<' '<<enums[m].size()<<'\n';
					for(size_t l=0; l<enums[m].size(); ++l) {
						int e = enums[m][l];
						int t = otherNode(e, m);
//						cout<<"other: "<<t<<' '<<s<<' '<<found[t]<<'\n';
						if (used[t]==s) continue;
						if (p>0) {
							cols[++z] = E*p + e;
							row[z] = 1;
						}
						cols[++z] = e;
						row[z] = 1;
					}
				}
#if 0
				cout<<" adding subtour cut "<<p<<' '<<cur.size()<<' '<<s<<' '<<z<<' '<<csum<<'\n';
				for(size_t i=0; i<cur.size(); ++i) cout<<cur[i]<<' ';cout<<'\n';
#endif
//				assert(z);
				if (!z) break;
				double viol = (p==0 ? 1 : 2) - csum;
				if (viol > maxViol) {
					maxViol = viol;
					if (rownum < 0) {
						rownum = glp_add_rows(lp, 1);
						++added;
					}
					glp_set_row_bnds(lp, rownum, GLP_LO, p==0 ? 1 : 2, 0);
					glp_set_mat_row(lp, rownum, z, &cols[0], &row[0]);
				}
			} else if (weakCuts && p==0 && s>1 && isum > (cur.size()-1)*csum + EPS) {
#if 1
				int z=0;
				for(size_t k=0; k<cur.size(); ++k) {
					int m = cur[k];
//					cout<<"conn from "<<m<<": "<<enums[m].size()<<'\n';
					for(size_t l=0; l<enums[m].size(); ++l) {
						int e = enums[m][l];
						int t = otherNode(e, m);
						if (used[t]==s && t<m) continue;
						cols[++z] = e;
						row[z] = used[t]==s ? -1 : (int)cur.size()-1;
//						cout<<"edge "<<m<<' '<<t<<' '<<vals[e]<<' '<<used[t]<<' '<<found[t]<<" ; "<<row[z]<<'\n';
					}
				}
//				cout<<"adding main path subtour cut "<<s<<' '<<cur.size()<<' '<<csum<<' '<<isum<<' '<<z<<' '<<cur<<'\n';

				double viol = isum - (cur.size()-1)*csum;
				if (viol > maxViol) {
					maxViol = viol;
					if (rownum < 0) {
						rownum = glp_add_rows(lp, 1);
						++added;
					}
					glp_set_row_bnds(lp, rownum, GLP_LO, 0, 0);
					glp_set_mat_row(lp, rownum, z, &cols[0], &row[0]);
				}
#endif
			}
		}
//		if (p==0 && s>1) cout<<"nodes explored: "<<cur.size()<<'\n';
	}
	if (added) cout<<"added "<<added<<" subtour cuts from "<<p<<'\n';
	return added;
}
int addSubtourConstraints()
{
	int added = 0;
	for(int i=0; i<=K; ++i) {
		added += addSubtoursSingle(i);
	}
	if (!added) added = addSubtoursSingle(0, 1);
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
			edists.push_back(edgeDist[i][j]);
		}
	}
	E = edges.size() - 1;
	lp = glp_create_prob();
	int vs = E*(K+1);
	glp_add_cols(lp, vs + robustOpt);
	for(int i=1; i<=vs; ++i) {
		int up = 2;
		if ((singleDir && i<=E) || singleDirAll) up=1;
		glp_set_col_bnds(lp, i, GLP_DB, 0, up);
	}

	cols.resize(vs+1);
	row.clear();
	row.resize(vs+1, 1);

	glp_set_obj_dir(lp, GLP_MIN);
	for(int i=1; i<=E; ++i) {
		glp_set_obj_coef(lp, i, edists[i] * LENGTH_FACTOR);
//		cout<<"lol "<<edges[i].first<<' '<<edges[i].second<<' '<<edists[i] * LENGTH_FACTOR<<'\n';
	}
	if (robustOpt) {
		glp_set_col_bnds(lp, 1+vs, GLP_LO, 0, 0);
		glp_set_obj_coef(lp, 1+vs, 1);
		int r = glp_add_rows(lp, K);
		for(int i=0; i<K; ++i, ++r) {
			int z=0;
			for(int j=1; j<=E; ++j) {
				cols[++z] = E*(1+i) + j;
//				row[z] = -edists[j] * probs[i];
				row[z] = -edists[j];
//				cout<<"k "<<i<<' '<<j<<' '<<cols[z]<<'\n';
			}
			cols[++z] = 1+vs;
			row[z] = 1;
//			cout<<"asd "<<i<<' '<<z<<'\n';
			glp_set_row_bnds(lp, r, GLP_LO, 0, 0);
			glp_set_mat_row(lp, r, z, &cols[0], &row[0]);
		}
	} else {
		for(int i=1; i<=K; ++i) {
			for(int j=1; j<=E; ++j) {
				glp_set_obj_coef(lp, E*i + j, edists[j] * probs[i-1]);
	//			cout<<"setting coeff "<<i<<' '<<j<<": "<<edists[j]*probs[i-1]<<' '<<edists[j]<<'\n';
			}
		}
	}

	{
		// Bounds for start and end nodes
		int r = glp_add_rows(lp, 2);
		for(int i=0; i<2; ++i, ++r) {
			for(size_t k=0; k<enums[i].size(); ++k) {
				cols[1+k] = enums[i][k];
				row[1+k] = 1;
//				cout<<"v "<<i<<' '<<k<<' '<<enums[i][k]<<'\n';
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
				row[z] = 1;
				cols[++z] = enums[t][k] + E*(1+i);
				row[z] = 1;
			}
			glp_set_row_bnds(lp, r, GLP_LO, 2, 0);
			glp_set_mat_row(lp, r, z, &cols[0], &row[0]);
		}
	}

//	glp_write_lp(lp, 0, "/dev/stdout");

	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.meth = GLP_PRIMAL;
	parm.presolve = 1;
	parm.msg_lev = GLP_MSG_ERR;

	vals.resize(1 + vs);
	do {
		int r = glp_simplex(lp, &parm);
		if (r || glp_get_status(lp)!=GLP_OPT) {
			cout<<"FAILED SOLVING ROUTE RELAXATION: "<<r<<' '<<glp_get_status(lp)<<'\n';
			dumpGraph(0);
			abort();
		}
		parm.meth = GLP_DUAL;
		parm.presolve = 0;

		cout<<"lower bound: "<<glp_get_obj_val(lp)<<'\n';

		for(int i=1; i<=vs; ++i)
			vals[i] = glp_get_col_prim(lp, i);
	} while(addConstraints());

	dumpGraph(0);

	double res = glp_get_obj_val(lp);

	glp_delete_prob(lp);

	return res;
}
