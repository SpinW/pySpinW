# Docstring Overview Improvements Needed: `tb2j.py`

## `TB2J_Input.from_exchange_out` (line 206), `TB2J_Input.from_text` (line 223), `TB2J_Input.from_pickle` (line 239)
- **Current overview**: `Read input from an exchange.out TB2J file.` / `Read input from text string containing the contents of an exchange.out TB2J file.` / `Read input from a TB2J.pickle file.`
- **Issue**: These classmethods appear intended as factory constructors returning a `TB2J_Input` instance, but they return the result of `cls()._parse_text(...)` or `cls()._parse_pickle(...)`, which mutate the created instance in place and return `None`.
- **Suggested clarification needed**: Clarify whether `_parse_text` / `_parse_pickle` (or the classmethods themselves) should return the constructed `TB2J_Input` instance.
