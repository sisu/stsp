#ifndef AREA_HPP
#define AREA_HPP

#include <istream>
#include <map>
#include <vector>
#include <string>
#include "Vector.hpp"
#include "Rect.hpp"

struct Area {
	double W,H;
	Vec2 startV, endV;
	std::vector<Rect> rects;
	std::map<std::string,Vec2> items;

	void read(std::istream& in);
	void readSimple(std::istream& in);
};

#endif
