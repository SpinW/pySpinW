# PySpinW will have a backwards compatibility layer

## Context

A clean class structure for pySpinW will entail making some relatively large changes to the intended use, this
might be off-putting to some users, and will break compatibility with existing scripts.

One big change is that there won't be a single big, stateful class, but lots of smaller stateless classes.

However, it would not be difficult to maintain existing syntax by putting all the new stuff into a stateful class that
closely resembles the original matlab.

## Decision

Create a class that mimics the original spinW interface, backed by the new objects. 

A subtle choice here would be to deliberately omit functionality from this class, and force "power users" to
migrate.

## Status

Proposed


## Consequences

* Will take time to do
* Might prevent users from migrating, and become permanent
