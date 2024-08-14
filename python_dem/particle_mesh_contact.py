import pygmsh

import meshio

with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        mesh_size=1.9,
    )
    mesh = geom.generate_mesh()

# mesh.points, mesh.cells, ...
mesh.write("out.vtk")
# mesh.write("out.obj")
# mesh.write("out.vtu")


mesh = meshio.read(
    "out.vtk"  # string, os.PathLike, or a buffer/open file
    # file_format="stl",  # optional if filename is a path; inferred from extension
    # see meshio-convert -h for all possible formats
)
