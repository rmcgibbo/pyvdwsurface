#ifndef _dotsphere_h
#define _dotsphere_h

#include <vector>
#include "Vec3.h"
using namespace OpenMM;
using namespace std;

vector<Vec3>
vdw_surface(vector<Vec3> coordinates, vector<string> elements,
            double scale_factor, double density);

#endif
                            