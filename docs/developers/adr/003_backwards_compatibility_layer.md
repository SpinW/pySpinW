# PySpinW will have a backwards compatibility layer

## Context

A clean class structure for pySpinW will entail making some relatively large changes to the intended use, this
might be off-putting to some users, and will break compatibility with existing scripts.

One big change is that there won't be a single big, stateful class, but lots of smaller stateless classes.

However, it would not be difficult to maintain existing syntax by putting all the new stuff into a stateful class that
closely resembles the original matlab.

## Decision

Create a class that mimics the original spinW interface, backed by the new objects. 
In particular, this class will only expose functionality which is currently in the Matlab version (and will omit new functionality written for Python) in order to encourage users to migrate.

## Status

Accepted

## Consequences

* Will take time to do
* Might prevent users from migrating, and become permanent
