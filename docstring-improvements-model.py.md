# Docstring Overview Improvements Needed: `model.py`

## `Model.render_back` (line 69)
- **Current overview**: `Render the back faces using wireframe - used for making selections.`
- **Issue**: The docstring states that back faces are rendered "using wireframe", whereas the implementation sets `glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)` with `glCullFace(GL_FRONT)` (filled back-face rendering for outline/selection silhouettes).
- **Suggested clarification needed**: Confirm whether the docstring should say "Render the back faces as filled triangles (front-face culled) for selection outlines" instead of "using wireframe".
