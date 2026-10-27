# Docstring Overview Improvements Needed: `anisotropy_shader.py`

## Module `pyspinw.gui.rendering.anisotropy_shader` (line 1)
- **Current overview**: `Shader for most of the objects in the viewer.`
- **Issue**: The overview appears to be copy-pasted from `object_shader.py` and does not describe anisotropy rendering or the Fresnel shader effect used here.
- **Suggested clarification needed**: Clarify that this module provides the shader program (`object_vertex` / `fresnel_fragment`) used to render magnetic anisotropy surfaces/ellipsoids in the 3D viewer.

## `AnisotropyShader` (line 10)
- **Current overview**: `Shader for most of the objects in the viewer.`
- **Issue**: Copy-pasted from `ObjectShader`; does not explain what `AnisotropyShader` specifically renders.
- **Suggested clarification needed**: Specify that `AnisotropyShader` renders anisotropy visualisations using a Fresnel fragment shader.
