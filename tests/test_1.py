import sys
# sys.path.append("/home/lo/projects/stVenant/build/src")

import meshpy.triangle as triangle
import stVenant
import numpy as np

circle_phi = np.linspace(0, 2 * np.pi, 100)[0:-1]

circle_x = np.cos(circle_phi) * 2
circle_y = np.sin(circle_phi)

points = np.array([circle_x, circle_y]).T.tolist()
bemCase = stVenant.TorsionBemCase([points])
bemCase.run()

print("number of elements: ", len(points))
print(bemCase.get_torsion_moment())
print(bemCase.get_polar_moment())

edges = np.array(range(len(points)))
edges = np.array([edges, edges + 1]).T
edges[-1, -1] = 0

info = triangle.MeshInfo()
info.set_points(points)
info.set_facets(edges)

mesh = triangle.build(info, max_volume=0.01, min_angle=25)
mesh_points = list(mesh.points)
triangles = list(mesh.elements)

femCase = stVenant.TorsionFemCase(mesh_points, triangles)
femCase.run()


print("number of elements: ", len(triangles))
print(femCase.get_torsion_moment())
print(femCase.get_polar_moment())

