pyvdwsurface
============

This package provides a single python function for computing points on the Van der Waals surface of a molecule.

![](https://raw.githubusercontent.com/rmcgibbo/pyvdwsurface/master/example/example.png)


```
def vdwsurface(coordinates, elements, scale_factor=1, density=1):
    Compute evenly-spaced points on the VDW surface of a molecule
    
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
        of surface area. 1.0 is the default recommended by Kollman & Singh
        for RESP calculations -- other uses may want higher densities.
    
    Returns
    -------
    points : np.ndarray, shape=(n_suface_points, 3)
        The cartesian coordinates of the surface points, in units of ANGSTROMS
```
