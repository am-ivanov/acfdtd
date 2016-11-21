#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <string>

enum Boundary {
	FREE_BOUNDARY,
	OPEN_BOUNDARY
};

struct Source {
	int x, y, z;
	std::vector<double> val;
};

struct Receiver {
	int x, y, z;
};

template <typename T>
inline T sqr(T x) {
	return x * x;
}

#endif // COMMON_H
