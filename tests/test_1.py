import meshpy.triangle as triangle
import stVenant
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import paraEigen as eigen


circle_phi = np.linspace(0, 2 * np.pi, 50)[0:-1]

circle_x = np.cos(circle_phi)
circle_y = np.sin(circle_phi) * 0.5

points = np.array([circle_x, circle_y]).T.tolist()
bemCase = stVenant.TorsionBemCase([points])
bemCase.run()

print("number of elements: ", len(points))
print("torsion-moment", bemCase.get_torsion_moment())
print("polar-moment", bemCase.get_polar_moment())

x = np.linspace(-1, 1, 100)
y = np.linspace(-0.6, 0.6, 100)
X, Y = np.meshgrid(x, y)
w = [[bemCase.get_w([xi, yi]) for xi in x] for yi in y]
w_min = min([wj for wi in w for wj in wi])
w_max = max([wj for wi in w for wj in wi])

plt.subplot(2, 1, 1)
plot = plt.contourf(x, y, w)
plt.colorbar(plot, orientation='vertical')

edges = np.array(range(len(points)))
edges = np.array([edges, edges + 1]).T
edges[-1, -1] = 0

info = triangle.MeshInfo()
info.set_points(points)
info.set_facets(edges)

mesh = triangle.build(info, max_volume=0.01)
mesh_points = list(mesh.points)
triangles = list(mesh.elements)

femCase = stVenant.TorsionFemCase(mesh_points, triangles)
femCase.run()


print("number of elements: ", len(triangles))
print("torsion-moment", femCase.get_torsion_moment())
print("polar-moment", femCase.get_polar_moment())

x, y = np.array(mesh_points).T
triang = mtri.Triangulation(x, y, triangles)

plt.subplot(2, 1, 2)
plt.triplot(triang, lw=0.5, color='0')
plot = plt.tricontourf(triang, femCase.get_w())
plt.colorbar(plot, orientation='vertical')

plt.show()
