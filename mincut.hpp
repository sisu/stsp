#ifndef MINCUT_HPP
#define MINCUT_HPP

#include <vector>
#include <utility>

double minCutFromTo(int s, int t, const std::vector<std::vector<int> >& enums, const std::vector<std::pair<int,int> >& edges, const std::vector<double>& ecaps, std::vector<int>& res);

#endif
