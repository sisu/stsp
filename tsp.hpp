#ifndef TSP_HPP
#define TSP_HPP

#include <vector>
std::vector<int> TSP(const std::vector<int>& p);
void initTSP();

extern std::vector<std::vector<std::vector<int> > > itemPath;
extern std::vector<std::vector<int> > startPath, endPath;

#endif
