#ifndef MINCUT_HPP
#define MINCUT_HPP

#include <vector>
#include <utility>

double minCutFromTo(int s, int t, const std::vector<std::vector<int> >& enums, const std::vector<std::pair<int,int> >& edges, const std::vector<double>& ecaps, std::vector<int>& res);

double minCutFrom(int start, int start2, const std::vector<int>& ends, std::vector<std::vector<std::pair<int,double> > >& edges, std::vector<int>& res);

#endif
