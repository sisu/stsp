#ifndef UTIL_HPP
#define UTIL_HPP

#include <ostream>
#include <vector>
#include <utility>
#include <tr1/functional>

template<class A, class B>
std::ostream& operator<<(std::ostream& s, std::pair<A,B> p)
{
	s<<'('<<p.first<<','<<p.second<<')';
	return s;
}
template<class T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v)
{
	s<<'[';
	for(size_t i=0; i<v.size(); ++i) {
		if (i) s<<',';
		s<<v[i];
	}
	s<<']';
	return s;
}

namespace std{namespace tr1{
template<class A,class B>
struct hash<pair<A,B> > : unary_function<pair<A,B>,size_t> {
	hash<A> a;
	hash<B> b;
	std::size_t operator()(const pair<A,B>& p) const {
		return a(p.first)*33331 + b(p.second);
	}
};
}}

#endif
