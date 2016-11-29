#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_FLOAT
typedef float real_t;
#else

typedef double real_t;
#endif

typedef long int_t;

enum Side { L = 0, R = 1, LR = 2 };

enum Boundary {
	FREE_BOUNDARY,
	OPEN_BOUNDARY
};

struct InterpPoint {
	real_t x, y, z;
	real_t c[LR][LR][LR]; // coefficient for bilinear interpolation
	int_t i, j, k; // coordinate of the left down bot node
};

struct Source : public InterpPoint {
	std::vector<real_t> val;
};

struct Receiver : public InterpPoint {
};

template <typename T>
inline T sqr(T x) {
	return x * x;
}

#endif // COMMON_H
