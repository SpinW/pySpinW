import numpy as np
from OpenGL.GL import *

class IntegerBuffer:
    def __init__(self):

        WIDTH, HEIGHT = 1024, 1024

        # Create integer texture
        id_tex = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, id_tex)
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_R32UI,
            WIDTH,
            HEIGHT,
            0,
            GL_RED_INTEGER,
            GL_UNSIGNED_INT,
            None
        )

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)

        # Create framebuffer
        self.fbo = glGenFramebuffers(1)
        glBindFramebuffer(GL_FRAMEBUFFER, self.fbo)

        glFramebufferTexture2D(
            GL_FRAMEBUFFER,
            GL_COLOR_ATTACHMENT0,
            GL_TEXTURE_2D,
            id_tex,
            0
        )

        # Depth buffer
        self.depth = glGenRenderbuffers(1)
        glBindRenderbuffer(GL_RENDERBUFFER, self.depth)
        glRenderbufferStorage(
            GL_RENDERBUFFER,
            GL_DEPTH_COMPONENT24,
            WIDTH,
            HEIGHT
        )

        glFramebufferRenderbuffer(
            GL_FRAMEBUFFER,
            GL_DEPTH_ATTACHMENT,
            GL_RENDERBUFFER,
            self.depth
        )

        assert glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE
        glBindFramebuffer(GL_FRAMEBUFFER, 0)

    def use(self, width: int, height: int):
        glBindFramebuffer(GL_FRAMEBUFFER, self.fbo)
        glViewport(0, 0, width, height)

        glDisable(GL_BLEND)  # REQUIRED
        glDisable(GL_DITHER)  # Recommended

        glClearBufferuiv(GL_COLOR, 0, [0])  # Clear integer buffer
        glClear(GL_DEPTH_BUFFER_BIT)

    def get_image(self, width, height):
        glBindFramebuffer(GL_FRAMEBUFFER, self.fbo)

        # Allocate numpy array
        img = np.empty((height, width), dtype=np.uint32)

        glReadPixels(
            0, 0,
            width, height,
            GL_RED_INTEGER,
            GL_UNSIGNED_INT,
            img
        )

        glBindFramebuffer(GL_FRAMEBUFFER, 0)

        # OpenGL origin is bottom-left
        return np.flipud(img)
