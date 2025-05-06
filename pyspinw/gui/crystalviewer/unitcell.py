from pyspinw.gui.crystalviewer.GL.models import WireModel


class UnitCellGraphics(WireModel):
    cube_vertices = [
        (-0.5, -0.5, -0.5),
        (-0.5, -0.5, 0.5),
        (-0.5, 0.5, -0.5),
        (-0.5, 0.5, 0.5),
        (0.5, -0.5, -0.5),
        (0.5, -0.5, 0.5),
        (0.5, 0.5, -0.5),
        (0.5, 0.5, 0.5)
    ]

    cube_edges = [
        (0, 1),  # Front face
        (1, 5),
        (5, 4),
        (4, 0),
        (2, 3),  # Back face
        (3, 7),
        (7, 6),
        (6, 2),
        (1, 3),  # between faces
        (5, 7),
        (4, 6),
        (0, 2)
    ]
    
    def __init__(self,
                 colors: Optional[ColorSpecification]=None,
                 edge_colors: Optional[ColorSpecification]=None):

        super().__init__(
            vertices=Cube.cube_vertices,
            edges=Cube.cube_edges,
            triangle_meshes=Cube.cube_triangles,
            edge_colors=edge_colors,
            colors=colors)

        self.vertices = Cube.cube_vertices
        self.edges = Cube.cube_edges

        if edge_colors is None:
            self.wireframe_render_enabled = False
            self.edge_colors = []
        else:
            self.wireframe_render_enabled = True
            self.edge_colors = edge_colors

        if colors is None:
            self.solid_render_enabled = False
            self.face_colors = []
        else:
            self.solid_render_enabled = True
            self.face_colors = colors
