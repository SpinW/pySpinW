# Docstring Overview Improvements Needed: `displayoptions.py`

## `DisplayOptions` (line 8)
- **Current overview**: `Options for how things should display.`
- **Issue**: The overview is very brief, and the field `show_atoms_not_spins` has a potentially confusing name compared to its usage in `CrystalViewerWidget.paintGL` (`show_arrow = site.is_magnetic and self.display_options.show_atoms_not_spins`).
- **Suggested clarification needed**: Clarify the intended semantics of `show_atoms_not_spins` (whether `True` means showing spin arrows or atomic spheres) and expand the class overview.
