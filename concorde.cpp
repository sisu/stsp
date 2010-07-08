#include "concorde.hpp"
#include <iostream>
using namespace std;

void ConcordeTSP::init() {
	dat = new CCdatagroup;
	CCutil_init_datagroup(dat);

	int N = dists.size();
	CCtsp_init_cutpool(&N, 0, &pool);

	rstate = new CCrandstate;
	CCutil_sprand(1, rstate);

	sel = new CCtsp_cutselect;
	CCtsp_init_cutselect(sel);

	int K = N*(N-1)/2;
	elist = new int[2*K];
	elens = new int[K];

	int k=0;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j) {
			elist[2*k] = j;
			elist[2*k+1] = i;
			elens[k] = 0;
			k++;
		}

	ptour = new int[N];
	for(int i=0; i<N; ++i) ptour[i] = i;

	CCutil_datagroup_perm(N, dat, ptour);
}

static char name[] = "tsp";
double ConcordeTSP::calc() {
	int N = dists.size();
	if (N < 3) return 0;
	int k=0;

	int add=100;
	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			dists[i][j] += add;
	dists[0][1] = 0;

	for(int i=0; i<N; ++i)
		for(int j=0; j<i; ++j)
			elens[k++] = dists[i][j];
	int K = N*(N-1)/2;
	bool silent = 1;
	int r = CCtsp_init_lp(&tsp, name, -1, NULL, N, dat, K, elist, elens, K, elist, elens, true, ptour, 1e100, pool, NULL, silent, rstate);
	if (r) cout<<"FAIL\n";
	CCtsp_cutselect_set_tols(sel, tsp, 1, silent);
//	CCtsp_cutting_loop(tsp, sel, 1, silent, rstate);
	CCtsp_subtour_loop(tsp, silent, rstate);

	return tsp->lowerbound - (N-1)*add;
}
