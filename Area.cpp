#include "Area.hpp"
#include <cstdlib>
using namespace std;

void Area::read(istream& in)
{
	in>>W>>H;
	in>>startV.x>>startV.y>>endV.x>>endV.y;

	string s;
	while(in>>s) {
		if (isdigit(s[0])) {
			Rect r;
			r.x1 = atof(s.c_str());
			in>>r.y1>>r.x2>>r.y2;
			rects.push_back(r);
		} else {
			double x,y;
			in>>x>>y;
			items[s] = Vec2(x,y);
		}
	}
}
void Area::readSimple(istream& in)
{
	in>>W>>H;
	in>>startV.x>>startV.y>>endV.x>>endV.y;

	int n;
	in>>n;
	rects.resize(n);
	for(int i=0; i<n; ++i) {
		Rect& r=rects[i];
		in>>r.x1>>r.y1>>r.x2>>r.y2;
	}
	int k;
	in>>k;
	string s;
	for(int i=0; i<k; ++i) {
		double x,y;
		in>>s>>x>>y;
		items[s] = Vec2(x,y);
	}
}
