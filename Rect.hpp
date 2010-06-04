#ifndef RECT_HPP
#define RECT_HPP

#include "Vector.hpp"

struct Rect {
	double x1,y1,x2,y2;

	int outcode(Vec2 v) {
		return (v.x<=x1 ? 1 : (v.x>=x2 ? 2 : 0))
			| (v.y<=y1 ? 4 : (v.y>=y2 ? 8 : 0));
	}
};

inline bool hitsRect(Vec2 a, Vec2 b, Rect r) {
	double x1=r.x1, y1=r.y1;
	double x2=r.x2, y2=r.y2;
	
	int o1=r.outcode(a);
	if (o1==0) return true;
	int o2=r.outcode(b);
	if (o2==0) return true;
	if ((o1&o2)!=0) return false;

	int ou = o1|o2;
	if ((ou&3)==0 || (ou&12)==0) return true;

	double d = (b.y-a.y) / (b.x-a.x);

	int ac=0, bc=0;
	if ((ou&1)!=0) {
		double ay = a.y + d*(x1-a.x);
		ac = ay<y1 ? 4 : (ay>y2 ? 8 : 0);
		if (ac==0) return true;
	} else ac |= ((o1&3)==0 ? o1 : o2) & 12;
	if ((ou&2)!=0) {
		double by = a.y + d*(x2-a.x);
		bc = by<y1 ? 4 : (by>y2 ? 8 : 0);
		if (bc==0) return true;
	} else bc |= ((o1&3)==0 ? o1 : o2) & 12;
	return (ac|bc)==12;
}

#endif
