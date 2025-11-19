Spacegroup Lookup
=================

Within the lookup for space groups there should be entries for:
 * The "international short" symbol
   * The name followed by the choice string in square brackets
   * The name followed by the choice string in round brackets 
 * The "international full" symbol


Unique axis choices
-------------------

The unique axis is specified by in long form

Numeric Choices
---------------

If there are symmetry choices, they can be specified by a colon then the choice.
An entry for a group with numeric choices with the choice omitted should default to choice 1.

Permutations
------------

For most groups, the permutations are implicit in the full name.

Some groups that have alternate settings specified by permutation have identical full names, specifically 
groups 67 and 68. These should have entries supplemented by the permutation in round brackets, and in square
brackets.

Rhombohedral and Hexagonal
--------------------------

Groups with full names that start with 'R' can have either a rhombohedral or hexagonal setting.
These should be suffixed with 'R' or 'H'
There should be an entry without an 'R' or 'H' suffix that defaults to the hexagonal setting

Other considerations
--------------------

Lookup should be case and space insensitive.
Tests for this not causing collisions should be performed.




Status
======

Pending