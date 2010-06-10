#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Rect.hpp"
using namespace std;

int H,W;

vector<string> area;

const double S=10;

int onum=1;
void addObs(double x, double y)
{
	if (x<=0 || y<=0 || x>=W*W || y>=S*H) return;
	cout<<'i'<<onum<<' '<<x<<' '<<y<<'\n';
	++onum;
}

int main(int argc, char* argv[])
{
	cin>>H>>W;
	string tmp;
	for(int i=0; i<H; ++i) cin>>tmp;
	area.resize(H);
	for(int i=0; i<H; ++i) cin>>area[i];

	cout<<S*W<<' '<<S*H<<'\n';

	int sx=-1,sy=-1,ex=-1,ey=-1;
	for(int i=0; i<H; ++i) {
		for(int j=0; j<W; ++j) {
			if (area[i][j]=='1') sx=j,sy=i;
			else if (area[i][j]=='2') ex=j,ey=i;
		}
	}
	cout<<S*sx<<' '<<S*sy<<'\n'<<S*ex<<' '<<S*ey<<'\n';

	vector<Rect> rs;
	for(int y=0; y<H; ++y) {
		for(int x=0; x<W; ++x) {
			if (area[y][x]=='4') {
				int w=0;
				while(x+w<W && area[y][x+w]=='4') ++w;
				bool ok=1;
				int h=0;
				do {
					++h;
					if (y+h>=H) break;
					for(int i=0; i<w; ++i) if (area[y+h][x+i]!='4') ok=0;
				} while(ok);
				cout<<S*x<<' '<<S*y<<' '<<S*(x+w)<<' '<<S*(y+h)<<'\n';
				Rect r={S*x,S*y,S*(x+w),S*(y+h)};
				rs.push_back(r);

				for(int i=0; i<h; ++i)
					for(int j=0; j<w; ++j)
						area[y+i][x+j]='0';
			}
		}
	}
	for(size_t i=0; i<rs.size(); ++i) {
		Rect r=rs[i];
		int w=r.x2-r.x1, h=r.y2-r.y1;
		int a = w/(2*S);
		int da = w/(a+1);
		for(int i=1; i<=a; ++i) {
			addObs(r.x1+da*i, r.y1-.5*S);
			addObs(r.x1+da*i, r.y2+.5*S);
		}
		int b = h/(2*S);
		int db = h/(b+1);
		for(int i=1; i<=b; ++i) {
			addObs(r.x1-.5*S, r.y1+db*i);
			addObs(r.x2+.5*S, r.y1+db*i);
		}
	}
}
