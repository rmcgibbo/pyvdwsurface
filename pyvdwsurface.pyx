import numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Vec3.h" namespace "OpenMM":
    cdef cppclass Vec3:
        Vec3()
        Vec3(double x, double y, double z)
        double operator[](int index) const

cdef extern from "vdwsurface.h":
    vector[Vec3] vdw_surface(vector[Vec3] coordinates, vector[string] elements,
                             double scale_factor, double density)


def vdwsurface(double[:, ::1] coordinates, elements, double scale_factor=1, double density=1):
    """vdwsurface(coordinates, elements, scale_factor=1, double density=1)
    
    Compute points on the VDW surface of a molecule
    
    Parameters
    ----------
    coordinates : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the nuclei, in units of ANGSTROMS
    elements : list, shape=(n_atoms)
        The element symbols (C, H, etc) for each of the nuceli
    scale_factor : float, default=1
        The points on the molecular surface are set at a distance of
        scale_factor * vdw_radius away from each of the atoms.
    density : float, default=1
        The (approximate) number of points to generate per square angstrom
        of surface area. 1.0 is the default recommended by Kollman & Singh.
    
    Returns
    -------
    points : np.ndarray, shape=(n_suface_points, 3)
        The cartesian coordinates of the surface points, in units of ANGSTROMS
    """
    cdef int i
    cdef vector[Vec3] coordinates_
    cdef vector[Vec3] surfpoints
    cdef double[::1] rvalue_
    assert coordinates.shape[1] == 3
    assert len(coordinates) == len(elements)
    
    for i in range(coordinates.shape[0]):
        coordinates_.push_back(
            Vec3(coordinates[i,0], coordinates[i,1], coordinates[i,2]))
        
    surfpoints = vdw_surface(coordinates_, elements, scale_factor, density)

    returnvalue = np.zeros((surfpoints.size(), 3))
    for i in range(surfpoints.size()):
        returnvalue[i,0] = surfpoints[i][0]
        returnvalue[i,1] = surfpoints[i][1]
        returnvalue[i,2] = surfpoints[i][2]
    return returnvalue
