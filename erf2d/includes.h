#ifndef __INCLUDES_H__
#define __INCLUDES_H__

#include <cmath>
#include <unordered_map>

#include "generic.h"
#include "unsteady_heat.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

typedef oomph::Vector<double> vecd;
typedef const vecd cvecd;
typedef const uint cuint;
typedef const double cdouble;

typedef void (*FluxFctPt)(const double &, const Vector<double>&, double &);

#endif