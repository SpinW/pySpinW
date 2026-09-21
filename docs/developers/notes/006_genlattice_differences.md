Differences between python `genlattice` and MATLAB
===============================================

We attempt to recreate `genlattice` as close as possible to the MATLAB version. 
However the different structure of pyspinw means that there will be some non-trivial differences.

Below are the MATLAB docs for `genlattice`, along with notes on the differences...



generates crystal lattice
=========================

### Syntax

`genlattice(obj,Name,Value)`

`R = genlattice(___)`

### Description

`genlattice(obj,Name,Value)` generates all necessary parameters to define
a lattice including space group symmetry and store the result it in the
[spinw.lattice] field.

`R = genlattice(___)` also returns the rotation matrix that
rotates the inpub basis vectors to the internal coordinate system.

Alternatively the lattice parameters can be given directly when the
[spinw] object is created using the `spinw(inpStr)` command, where struct
contains the fields with initial parameters, e.g.:
```
inpStr.lattice.lat_const = [3 3 4];
```

### Example

```
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym','P 6')
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym',168)
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym','-y,x-y,z; -x,-y,z','label','R -3 m')
```

The three lines are equivalent, both will create hexagonal lattice, with
$P6$ space group.

### Input

`obj`
: [spinw] object.

### Options

`angled`
: `[\\alpha, \\beta, \\gamma]` angles in \\deg, dimensions are $[1\times 3]$.

`angle`
: `[\\alpha, \\beta, \\gamma]` angles in radian, dimensions are $[1\times 3]$.

`lat_const`
: `[a, b, c]` lattice parameters in units defined in [spinw.unit] (with \\ang
  being the default), dimensions are $[1\times 3]$.

`spgr` or 'sym'
: Defines the space group. Can have the following values:

  * **space group label** string, name of the space group, can be any
    label defined in the [symmetry.dat] file.
  * **space group index** line number in the [symmetry.dat] file.
  * **space group operators** matrix with dimensions 
    $[3\times 4\times n_{op}]$.
  
  The [symmetry.dat] file stores definition of the 230 space groups in
  standard settings as it is in the [International Tables of Crystallography](http://it.iucr.org/A/).
  Additional lines can be added to the [symmetry.dat] file using the
  [swsym.add] function which later can be used in the `spgr` option.

  If the `spgr` option is 0, no symmetry will be used. The
  [spinw.gencoupling] function will determine the equivalent bonds based on
  bond length.
  
  Can also provide spacegroup and label (see below) in a cell e.g.
  {'-x,y,-z', 'P 2'}

`label`
: Optional label for the space group if the generators are given in the
  `spgr` option.

`bv`
: Basis vectors given in a matrix with dimensions of $[3\times 3]$, where
  each column defines a basis vector.

`origin`
: Origin for the space group operators, default value is `[0 0 0]`.

`perm`
: Permutation of the abc axes of the space group operators.

`nformula`
: Gives the number of formula units in the unit cell. It is used
  to normalize cross section in absolute units. Default value is 0, when
  cross section is normalized per unit cell.

### Output

`R`
: Rotation matrix that brings the input basis vector to the SpinW
  compatible form:
  ```
  BVspinw = R*BV
  ```

The result of the `spinw.genlattice` function is that `obj.lattice` field
will be changed based on the input, the lattice parameters are stored
directly and the optional space group string is converted into space
group operator matrices.

### See also

[spinw], [swsym.add], [swsym.operator], [spinw.gencoupling]
