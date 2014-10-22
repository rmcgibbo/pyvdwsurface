import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from pyvdwsurface import vdwsurface

atoms = np.array([[-1,-1,0,], [-1,1,0], [1,-1,0], [1,1,0],], dtype=float)
points = vdwsurface(atoms, elements=['C']*len(atoms), density=20)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0], points[:,1], points[:,2], marker='o')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.set_zlim(-2,2)
plt.savefig('example.png')
