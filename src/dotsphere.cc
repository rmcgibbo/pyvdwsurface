//  This code is adapted from GROMACS
//
//  Copyright (c) 1991-2000, University of Groningen, The Netherlands.
//  Copyright (c) 2001-2007, The GROMACS development team,
//  check out http://www.gromacs.org for more information.
//  Copyright (c) 2012,2013, by the GROMACS development team, led by
//  David van der Spoel, Berk Hess, Erik Lindahl, and including many
//  others, as listed in the AUTHORS file in the top-level source
//  directory and at http://www.gromacs.org.
//
//  GROMACS is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public License
//  as published by the Free Software Foundation; either version 2.1
//  of the License, or (at your option) any later version.
//
//  GROMACS is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////////
// Two routines for computing a set of dots (approximately) unfiromly
// distributed on the unit sphere by repeated truncation of icosahedrons.
//
// The two methods, dotsphere1 and dotsphere2, use slightly different procedures
// which causes them to be capable of yielding different numbers of points.
//////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <set>

#include "Vec3.h"
#include "dotsphere.h"

#define TORAD(A)     ((A)*0.017453293)
#define DP_TOL     0.001

// Number of steps of electrostatic refinement
#define N_REFINE_STEPS 25

using namespace std;
using namespace OpenMM;

static double R_H = sqrt(1.0 - 2.0*cos(TORAD(72.)))/(1.-cos(TORAD(72.)));

Vec3 divarc(const Vec3& xyz1, const Vec3& xyz2, int div1, int div2) {
    // Divide an arc based on the great circle
    double xd = xyz1[1]*xyz2[2] - xyz2[1]*xyz1[2];
    double yd = xyz1[2]*xyz2[0] - xyz2[2]*xyz1[0];
    double zd = xyz1[0]*xyz2[1] - xyz2[0]*xyz1[1];
    double dd = sqrt(xd*xd + yd*yd + zd*zd);
    if (dd < DP_TOL)
        throw "_divarc: rotation axis of length";

    double d1 = xyz1.norm();
    if (d1 < sqrt(0.5))
        throw "_divarc: vector 1 of sq.length too small";

    double d2 = xyz2.norm();
    if (d2 < sqrt(0.5))
        throw "_divarc: vector 2 of sq.length too small";

    double phi = sin(dd / sqrt(d1*d2));
    phi = phi * static_cast<double>(div1) / static_cast<double>(div2);
    double sphi = sin(phi);
    double cphi = cos(phi);
    double s = (xyz1[0]*xd + xyz1[0]*yd + xyz1[2]*zd) / dd;

    double x = xd*s*(1.0-cphi)/dd + xyz1[0] * cphi + (yd*xyz1[2] - xyz1[1]*zd)*sphi/dd;
    double y = yd*s*(1.0-cphi)/dd + xyz1[1] * cphi + (zd*xyz1[0] - xyz1[2]*xd)*sphi/dd;
    double z = zd*s*(1.0-cphi)/dd + xyz1[2] * cphi + (xd*xyz1[1] - xyz1[0]*yd)*sphi/dd;
    dd = sqrt(x*x + y*y + z*z);

    return Vec3(x/dd, y/dd, z/dd);
}

vector<Vec3> icosahedron_vertices(void) {
    // Compute the vertices of an icosahedron
    //
    // This code was adapted from GROMACS's nsc.c, distributed under the GNU
    // LGPL. See this file's header for the copyright information.
    //
    // Returns
    // -------
    // verts : np.ndarray, shape=(12, 3)
    //     Cartesian coordinates of the 12 verticles of a unit icosahedron.
    //

    double rg = cos(TORAD(72.))/(1.-cos(TORAD(72.)));

    std::vector<Vec3> verts;
    verts.push_back(Vec3(0.0, 0.0, 1.0));
    verts.push_back(Vec3(R_H*cos(TORAD(72.)), R_H*sin(TORAD(72.)), rg));
    verts.push_back(Vec3(R_H*cos(TORAD(144.)), R_H*sin(TORAD(144.)), rg));
    verts.push_back(Vec3(R_H*cos(TORAD(216.)), R_H*sin(TORAD(216.)), rg));
    verts.push_back(Vec3(R_H*cos(TORAD(288.)), R_H*sin(TORAD(288.)), rg));
    verts.push_back(Vec3(R_H, 0, rg));
    verts.push_back(Vec3(R_H*cos(TORAD(36.)), R_H*sin(TORAD(36.)), -rg));
    verts.push_back(Vec3(R_H*cos(TORAD(108.)), R_H*sin(TORAD(108.)), -rg));
    verts.push_back(Vec3(-R_H, 0, -rg));
    verts.push_back(Vec3(R_H*cos(TORAD(252.)), R_H*sin(TORAD(252.)), -rg));
    verts.push_back(Vec3(R_H*cos(TORAD(324.)), R_H*sin(TORAD(324.)), -rg));
    verts.push_back(Vec3(0, 0, -1));

    return verts;
}


vector<Vec3> dotsphere1(int density) {
    // Create a dot distribution over the unit shpere based on repeated
    // splitting and refining the arcs of an icosahedron.
    //
    // In general, by my visual inspection (RTM, August 2013), the dot
    // distributions produced by this method look worse than those produced by
    // the alternative procedure, `dotsphere_icos2`. This method generally
    // produces fewer points.
    //
    // Parameters
    // ----------
    // density : int
    //     Required number of dots on the unit sphere
    //
    // Returns
    // -------
    // dots : np.ndarray, shape=(N, 3), dtype=np.double
    //     Dots on the surface of the unit sphere. The number of dots will be
    //     at minimum equal to the `density` argument, but will be roughly two
    //     times larger.
    //
    // Notes
    // -----
    // This code was adapted from the function 'ico_dot_arc' in GROMACS's nsc.c,
    // distributed under the GNU LGPL. See this file's header for the copyright
    // information.
    //
    // See Also
    // --------
    // dotsphere_icos2 : acomplished the same goal, but based on splitting
    //     the faces. The two procedures are capable of yielding different
    //     number of points because of the different algorithms used.

    // calculate tessalation level
    double d;
    double a = sqrt((density - 2.0) / 10.0);
    int tess = static_cast<int>(ceil(a));

    vector<Vec3> vertices = icosahedron_vertices();

    if (tess > 1) {
        a = R_H*R_H*2.0*(1.0 - cos(TORAD(72.0)));
        // Calculate tessalation of icosahedron edges
        for (int i = 0; i < 11; i++) {
            for (int j = i+1; j < 12; j++) {
                d = (vertices[i] - vertices[j]).norm();
                if (fabs(a-d*d) > DP_TOL)
                    continue;
                for (int tl = 0; tl < tess; tl++)
                    vertices.push_back(divarc(vertices[i], vertices[j], tl, tess));
            }
        }
    }

    // Calculate tessalation of icosahedron faces
    for (int i = 0; i < 10; i++) {
        for (int j = i+1; j < 11; j++) {
            d = (vertices[i] - vertices[j]).norm();
            if (fabs(a-d*d) > DP_TOL)
                continue;

            for (int k = j+1; k < 12; k++) {
                double d_ik = (vertices[i] - vertices[k]).norm();
                double d_jk = (vertices[j] - vertices[k]).norm();
                if ((fabs(a - d_ik*d_ik) > DP_TOL) || (fabs(a - d_jk*d_jk) > DP_TOL))
                    continue;
                for (int tl = 1; tl < tess-1; tl++) {
                    Vec3 ji = divarc(vertices[j], vertices[i], tl, tess);
                    Vec3 ki = divarc(vertices[k], vertices[i], tl, tess);

                    for (int tl2 = 1; tl2 < tess-tl; tl2++) {
                        Vec3 ij = divarc(vertices[i], vertices[j], tl2, tess);
                        Vec3 kj = divarc(vertices[k], vertices[j], tl2, tess);
                        Vec3 ik = divarc(vertices[i], vertices[k], tess-tl-tl2, tess);
                        Vec3 jk = divarc(vertices[j], vertices[k], tess-tl-tl2, tess);

                        Vec3 xyz1 = divarc(ki, ji, tl2, tess-tl);
                        Vec3 xyz2 = divarc(kj, ij, tl, tess-tl2);
                        Vec3 xyz3 = divarc(jk, ik, tl, tl+tl2);

                        Vec3 x = xyz1 + xyz2 + xyz3;
                        vertices.push_back(x / x.norm());
                    }
                }
            }
        }
    }
    return vertices;
}

vector<Vec3> dotsphere2(int density) {
    // Create a dot distribution over the unit shpere based on repeated
    // truncating and refining the faces of an icosahedron.
    //
    // In general, by my visual inspection (RTM, August 2013), the dot
    // distributions produced by this method look "better" than those produced by
    // the alternative procedure, `dotsphere_icos1`. But this method tends to
    // produce more points.
    //
    // Parameters
    // ----------
    // density : int
    //     Required number of dots on the unit sphere
    //
    // Returns
    // -------
    // dots : np.ndarray, shape=(N, 3), dtype=np.double
    //     Dots on the surface of the unit sphere. The number of dots will be
    //     at minimum equal to the `density` argument, but will be roughly two
    //     times larger.
    //
    // Notes
    // -----
    // This code was adapted from the function 'ico_dot_dod' in GROMACS's nsc.c,
    // distributed under the GNU LGPL. See this file's header for the copyright
    // information.
    //
    // See Also
    // --------
    // dotsphere_icos1 : acomplished the same goal, but based on splitting
    //     the edges. The two procedures are capable of yielding different
    //     number of points because of the different algorithms used.

    double a = sqrt((density - 2.0) / 30.0);
    int tess = max(static_cast<int>(ceil(a)), 1);
    vector<Vec3> vertices = icosahedron_vertices();

    a = R_H*R_H * 2.0 * (1.0 - cos(TORAD(72.0)));
    
    // Dodecaeder vertices
    for (int i = 0; i < 10; i++) {
        for (int j = i+1; j < 11; j++) {
            double d = (vertices[i] - vertices[j]).norm();
            if (fabs(a-d*d) > DP_TOL)
                continue;
            for (int k = j+1; k < 12; k++) {
                double d_ik = (vertices[i] - vertices[k]).norm();
                double d_jk = (vertices[j] - vertices[k]).norm();
                if ((fabs(a - d_ik*d_ik) > DP_TOL) || (fabs(a - d_jk*d_jk) > DP_TOL))
                    continue;

                Vec3 x = (vertices[i] + vertices[j] + vertices[k]);
                vertices.push_back(x / x.norm());
            }
        }
    }

    if (tess > 1) {
        // square of the edge of an dodecaeder
        double adod = 4.0 * (cos(TORAD(108.)) - cos(TORAD(120.))) / (1.0 - cos(TORAD(120.)));
        // square of the distance of two adjacent vertices of ico- and dodecaeder
        double ai_d = 2.0 * (1.0 - sqrt(1.0 - a/3.0));
        // calculate tessalation of mixed edges
        for (int i = 0; i < 31; i++) {
            int j1 = 12;
            int j2 = 32;
            a = ai_d;
            if (i > 12) {
                j1 = i+1;
                a = adod;
            }
            for (int j = j1; j < j2; j++) {
                double d = (vertices[i] - vertices[j]).norm();
                if (fabs(a-d*d) > DP_TOL)
                    continue;
                for (int tl = 1; tl < tess; tl++) {
                    vertices.push_back(divarc(vertices[i], vertices[j], tl, tess));
                }
            }
        }
        
        // calculate tessalation of pentakisdodecahedron faces
        for (int i = 0; i < 12; i++) {
            for (int j = 12; j < 31; j++) {
                double d = (vertices[i] - vertices[j]).norm();
                if (fabs(ai_d-d*d) > DP_TOL)
                    continue;

                for (int k = j+1; k < 32; k++) {
                    double d_ik = (vertices[i] - vertices[k]).norm();
                    double d_jk = (vertices[j] - vertices[k]).norm();
                    if ((fabs(ai_d - d_ik*d_ik) > DP_TOL) || (fabs(adod - d_jk*d_jk) > DP_TOL))
                        continue;
                    for (int tl = 1; tl < tess-1; tl++) {
                        Vec3 ji = divarc(vertices[j], vertices[i], tl, tess);
                        Vec3 ki = divarc(vertices[k], vertices[i], tl, tess);
                        for (int tl2 = 1; tl2 < tess-tl; tl2++) {
                            Vec3 ij = divarc(vertices[i], vertices[j], tl2, tess);
                            Vec3 kj = divarc(vertices[k], vertices[j], tl2, tess);
                            Vec3 ik = divarc(vertices[i], vertices[k], tess-tl-tl2, tess);
                            Vec3 jk = divarc(vertices[j], vertices[k], tess-tl-tl2, tess);

                            Vec3 xyz1 = divarc(ki, ji, tl2, tess-tl);
                            Vec3 xyz2 = divarc(kj, ij, tl, tess-tl2);
                            Vec3 xyz3 = divarc(jk, ik, tl, tl+tl2);

                            Vec3 x = (xyz1 + xyz2 + xyz3);
                            vertices.push_back(x / x.norm());
                        }
                    }
                }
            }
        }
    }
    return vertices;
}

/****************************************************************************/
// Code for refining a dot distribution based on electrostatic repulsion
/****************************************************************************/

double get_coulomb_energy(const vector<Vec3>& points) {
    // Calculate the coulomb energy between a set of ponts
    double e = 0;
    for (size_t i = 0; i < points.size(); i++) {
        for (size_t j = i+1; j < points.size(); j++) {
            e += 1 / (points[i] - points[j]).norm();
        }
    }
    return e;
}

void get_coulomb_forces(vector<Vec3>& forces, const vector<Vec3>& points) {
  if (forces.size() != points.size())
      throw "Size Mismatch";
  for (size_t i = 0; i < forces.size(); i++) {
      forces[i][0] = 0;
      forces[i][1] = 0;
      forces[i][2] = 0;
  }

  for (size_t i = 0; i < points.size(); i++) {
      for (size_t j = i+1; j < points.size(); j++) {
          Vec3 r = points[i] - points[j];
          double l = r.norm();
          if (l < 1e-8)
              continue;
          Vec3 ff = r / (l*l*l);
          forces[i] += ff;
          forces[j] -= ff;
    }
  }
}

void refine_dotsphere(vector<Vec3>& vertices) {
    // Refine a dot distribution over the unit sphere using an iterative
    // electrostatic repulsion-type approach.

    // Parameters
    // ----------
    // vertices : vector<Vec3>
        // The initial vertices

    // Returns
    // -------
    // None. The vertices will be modified inplace
    double step = 0.005;
    if (vertices.size() > 100)
        step /= 50;

    vector<Vec3> forces(vertices.size());

    double e0 = get_coulomb_energy(vertices);
    for (size_t k = 0; k < N_REFINE_STEPS; k++) {
        get_coulomb_forces(forces, vertices);
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] + forces[i]*step;
            vertices[i] = vertices[i] / vertices[i].norm();
        }
        double e = get_coulomb_energy(vertices);

        if (e0 < e)
            step /= 2;
        e0 = e;

        if (step < 1e-8)
            break;
    }
}

vector<Vec3> dotsphere(int density) {
    // Create a dot distribution over the unit sphere, choosing the most
    // appropriate implementation based on the number of dots you request.
    //
    // Parameters
    // ----------
    // density : int
    //     Required number of dots on the unit sphere
    //
    
    unsigned int seed;
    int i1 = 1;
    int i2 = 1;
    vector<Vec3> vertices;

    while (10*i1*i1+2 < density)
        i1 += 1;

    while (30*i2*i2+2 < density)
        i2 += 1;
    
    // Use one of the two algorithms to generate the initial dots
    // they will give slightly too many.
    if (10*i1*i1-2 < 30*i2*i2-2) {
        vertices = dotsphere1(density);
    } else {
        vertices = dotsphere2(density);
    }

    if (static_cast<size_t>(density) < vertices.size()) {
        // Now lets throw out some of them
        set<size_t> keep;
        while (keep.size() < (size_t) density) {
            size_t v = rand_r(&seed) % vertices.size();
            keep.insert(v);
        }

        vector<Vec3> new_vertices;
        for (set<size_t>::iterator iter = keep.begin(); iter != keep.end(); iter++) {
            new_vertices.push_back(vertices[*iter]);
        }
        vertices = new_vertices;
    }
    refine_dotsphere(vertices);
    return vertices;
}
