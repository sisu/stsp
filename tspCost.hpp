#ifndef TSPCOST_HPP
#define TSPCOST_HPP

#include <vector>
double expectedTotalCost(const std::vector<int>& path);
double exactTotalCost(const std::vector<int>& path);
void initTSPCost(const std::vector<std::vector<int> >& ss, const std::vector<double> ps);

#endif
