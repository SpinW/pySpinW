import numpy as np
from OpenGL.GL import *

def read_int_fbo(fbo, width, height):
    glBindFramebuffer(GL_FRAMEBUFFER, fbo)

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