#pragma once

// useless :(
#include <stdfloat>

#ifndef A_TOL
#define A_TOL 1e-6
#endif

#ifndef R_TOL
#define R_TOL 0.0
#endif

#ifndef FPOINT_T
using fpoint_t = double;
#else
using fpoint_t = FPOINT_T;
#endif