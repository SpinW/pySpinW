# Docstring Overview Improvements Needed: `viewer.py`

## `_generate_unique_id` (line 43)
- **Current overview**: `Generate a unique ID for each site currently loaded.`
- **Issue**: The summary states that it generates an ID for each site, but in `Viewer.__init__` it is actually used to assign a unique ID to each `Viewer` window instance stored in `_VIEWERS`.
- **Suggested clarification needed**: Update the summary to state that it generates a unique ID for each active `Viewer` instance.
