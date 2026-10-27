# Docstring Overview Improvements Needed: `edge_cases.py`

## `add_extra_edge_lines` (line 44)
- **Current overview**: `When a line lives on the x/y/z=0 face/edges of the unit cube, make extra ones.` (plus bulleted extended summary)
- **Issue**: The extended summary contains typos and a factual error in the coordinate list (`an not an edge` instead of `and not an edge`, `2 two` instead of `two`, and `[0,0], [0,1], [1,0] and [0,0]` instead of `[0,0], [0,1], [1,0] and [1,1]`).
- **Suggested clarification needed**: Fix the typos in the extended summary bullets so the four 2D hypercube edge replacements are listed as `[0,0], [0,1], [1,0], and [1,1]`.
