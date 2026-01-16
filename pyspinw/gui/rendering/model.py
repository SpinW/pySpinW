class Model:

    def __init__(self):
        self._vertices_and_normals = None
        self._vao = None
        self._vbo = None


    @property
    def vao(self):
        return self._vao


    @property
    def vbo(self):
        return self._vbo


    def render(self, shader_program, model_matrix: np.ndarray, camera: Camera):
        pass