## pySpinW Classes with Fitable Parameters will Provide References to them through a property called `.parameters`

# Context

The fitting in spinW is reasonably limited. To do powerful and general fitting, we need a general way of specifying which 
parameters to fit.

# Decision

Objects in `parameters` will have:
* A name
* A public method for reading the parameter
* A private method for updating the parameter (though an immutable framework)
* Information about any fitting constraints that the parameter might require (e.g. positiveness, other bounds) 
* Private information that allows for checking for over-specification (see example below)

For example, the `parameters` for DM or diagonal coupling might be
* magnitude
* x_component
* y_component
* z_component

In a case like this, spinW will need to be aware that changing the magnitude and, say, x, will result in conflicts.

Fitting calls will work roughly like in the following sketch for a two parameter fitting:
```python
dm_coupling = DMCoupling(A, B, [1,1,1])
diag_coupling = DiagonalCoupling(A, B, [1,2,3])

...

instrument = Instrument.from_data(...)

...

instrument.fit(dm_coupling.parameters.magnitude, diag_coupling.parameters.magnitude, ...)


```



# Status

Proposed

# Advantages

PySpinW Needs a way of specifying fits that is relatively user-friendly and quick, but it also needs to be flexible.
This should achieve this compromise.