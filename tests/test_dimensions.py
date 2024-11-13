""" Tests for the check_dimensions decorator

Soooo many decorators, it must be Christmas.

It's November.

So yes.
"""

# Disable linting for unused arguments

# pylint: disable=W0613

import numpy as np
import pytest

from pyspinw.dimensionality import dimensionality_check, DimensionalityError

# Some lengths to check
vector_check_lengths = [0, 1, 3, 7]

#
# Groups of test cases: bad specifications, vector, matrix, vector-vector, crazy multitensor
#

# Note: I think I broke python with this one:

@pytest.mark.parametrize("bad_value", [3, 3.14, "a string", {"dict": "ionary"}])
def test_failed_definition_type(bad_value):
    """ Check silly inputs the decorator throw"""
    with pytest.raises(TypeError):
        @dimensionality_check(x=bad_value)
        def fun(x):
            pass


#
# Single parameter: vector cases
#

@dimensionality_check(x=(7,))
def fixed_length_vector(x):
    """ Test function that should be supplied with a length 7 array"""

@dimensionality_check(x=('n',))
def string_arbitrary_length_vector(x):
    """ Any vector should pass """

@dimensionality_check(x=(-1,))
def int_arbitrary_length_vector(x):
    """ Any vector should pass"""

@pytest.mark.parametrize("length, success", zip(vector_check_lengths, [False, False, False, True]))
def test_fixed_vector_parameter(length, success):
    """ Check fixed length vectors work """
    test_input = np.zeros((length,))

    # Check vector with right name
    if not success:
        with pytest.raises(DimensionalityError):
            fixed_length_vector(test_input)

        # Should work if keyword arguments are supplied instead
        with pytest.raises(DimensionalityError):
            fixed_length_vector(x=test_input)

    else:
        # TODO: check this is how to test a function just runs without error
        fixed_length_vector(test_input)
        fixed_length_vector(x=test_input)

    # The following this should throw
    #  * wrong variable name specified
    #  * wrong dimensionality

    with pytest.raises(Exception):

        # pylint: disable=no-value-for-parameter
        fixed_length_vector(y=test_input)
        # pylint: enable=no-value-for-parameter

    with pytest.raises(DimensionalityError):
        fixed_length_vector(np.zeros((4,5)))


@pytest.mark.parametrize("function", [string_arbitrary_length_vector,
                                      int_arbitrary_length_vector])
@pytest.mark.parametrize("length", vector_check_lengths)
def test_vector_parameter(function, length):
    """ Test arbitrary length vectors work, but still have to be vectors """
    test_input = np.zeros((length,))

    # Should just not fail
    function(test_input)
    function(x=test_input)

    # The following this should throw
    #  * wrong variable name specified
    #  * wrong dimensionality

    with pytest.raises(Exception):
        function(y=test_input)

    with pytest.raises(DimensionalityError):
        function(np.zeros((4,5)))


#
# Single parameter - matrix tests
#


@dimensionality_check(x=('n','n'))
def square_matrix(x):
    """ Constrained by having both dimensions of x equal to the same thing"""

@pytest.mark.parametrize('a', vector_check_lengths)
@pytest.mark.parametrize('b', vector_check_lengths)
def test_square_matrix_check(a, b):
    """ Check all square matrices pass, others don't"""
    mat = np.zeros((a, b))
    if a == b:
        square_matrix(mat)
    else:
        with pytest.raises(DimensionalityError):
            square_matrix(mat)

@dimensionality_check(x=('n','m'))
def any_matrix_all_names(x):
    """ Named but not constrained """

@dimensionality_check(x=('n',-1))
def any_matrix_name_number(x):
    """ half named but not constrained """

@dimensionality_check(x=(-1,'m'))
def any_matrix_number_name(x):
    """ half named but not constrained """

@dimensionality_check(x=(-1, -1))
def any_matrix_all_numbers(x):
    """ half named but not constrained """

@pytest.mark.parametrize('function', [ any_matrix_all_names,
                                       any_matrix_name_number,
                                       any_matrix_number_name,
                                       any_matrix_all_numbers])
@pytest.mark.parametrize('a', vector_check_lengths)
@pytest.mark.parametrize('b', vector_check_lengths)
def test_any_matrix(function, a, b):
    """ Unconstrained"""
    mat = np.zeros((a, b))
    function(mat)

@pytest.mark.parametrize('n', [1,2,3,4])
def test_any_matrix_dimensionality(n):
    """ Check that only 2D is allowed"""

    m = np.zeros(tuple([3]*n)) # 3, 3x3, 3x3x3 ...

    if n == 2:
        any_matrix_all_numbers(m)
    else:
        with pytest.raises(DimensionalityError):
            any_matrix_all_numbers(m)
@dimensionality_check(x=(1,7))
def one_by_seven_matrix(x):
    """ Input requires 1x7 matrix"""

@pytest.mark.parametrize('a', vector_check_lengths)
@pytest.mark.parametrize('b', vector_check_lengths)
def test_fixed_size_matrix(a, b):
    """ Check for a matrix with two fixed values"""
    mat = np.zeros((a, b))

    if a == 1 and b == 7:
        one_by_seven_matrix(mat)
    else:
        with pytest.raises(DimensionalityError):
            one_by_seven_matrix(mat)

@dimensionality_check(x=(7,-1))
def first_seven_matrix(x):
    """ First dimension is 7 long"""

@dimensionality_check(x=(-1, 7))
def second_seven_matrix(x):
    """ Second dimension is 7 long"""

@pytest.mark.parametrize('a', vector_check_lengths)
@pytest.mark.parametrize('b', vector_check_lengths)
def test_half_fixed(a, b):
    """ Checks for n-by-fixed and fixed-by-n"""
    mat = np.zeros((a, b))

    if a == 7:
        first_seven_matrix(mat)
    else:
        with pytest.raises(DimensionalityError):
            first_seven_matrix(mat)

    if b == 7:
        second_seven_matrix(mat)
    else:
        with pytest.raises(DimensionalityError):
            second_seven_matrix(mat)


#
# Vector-vector tests, only really need to check named combinations, and ommited variable
#

@dimensionality_check(x=('n',), y=('n', ))
def matched_vectors(x, y):
    """ Two vector inputs that must be the same length"""


@pytest.mark.parametrize('a', vector_check_lengths)
@pytest.mark.parametrize('b', vector_check_lengths)
def test_matched(a, b):
    """ Should give an error if called with different length vectors, run otherwise"""
    x = np.zeros((a, ))
    y = np.zeros((b, ))

    if a == b:
        matched_vectors(x, y)
    else:
        with pytest.raises(DimensionalityError):
            matched_vectors(x, y)

#
# combinations of multiple tensors
#
@dimensionality_check(x=('i','k'), y=('j', 'k'), z=('j', 'i'))
def triple_constraint(x, y, z):
    """ Three mutually constrained variables"""

only_three_lengths = [1, 3, 7]
@pytest.mark.parametrize('xa', only_three_lengths)
@pytest.mark.parametrize('xb', only_three_lengths)
@pytest.mark.parametrize('ya', only_three_lengths)
@pytest.mark.parametrize('yb', only_three_lengths)
@pytest.mark.parametrize('za', only_three_lengths)
@pytest.mark.parametrize('zb', only_three_lengths)
def test_triple(xa, xb, ya, yb, za, zb):
    """ Test a triple matrix input specification with interrelated constraints"""
    should_pass = xa == zb and xb == yb and ya == za

    x = np.zeros((xa, xb))
    y = np.zeros((ya, yb))
    z = np.zeros((za, zb))

    if should_pass:
        triple_constraint(x,y,z)
    else:
        with pytest.raises(DimensionalityError):
            triple_constraint(x,y,z)

@dimensionality_check(x=(3,'i','i','j'), y=('j', 7))
def tensor_1(x, y):
    """ A big tensor"""

@pytest.mark.parametrize('xa', only_three_lengths)
@pytest.mark.parametrize('xb', only_three_lengths)
@pytest.mark.parametrize('xc', only_three_lengths)
@pytest.mark.parametrize('xd', only_three_lengths)
@pytest.mark.parametrize('ya', only_three_lengths)
@pytest.mark.parametrize('yb', only_three_lengths)
def test_tensor_1(xa, xb, xc, xd, ya, yb):
    """ Test a function with big tensor, some self-consistency constraints, some interrelated, some fixed"""
    should_pass = xa == 3 and xb == xc and xd == ya and yb == 7

    x = np.zeros((xa, xb, xc, xd))
    y = np.zeros((ya, yb))

    if should_pass:
        tensor_1(x, y)
    else:
        with pytest.raises(DimensionalityError):
            tensor_1(x, y)

@dimensionality_check(x=('n','n','n','m'))
def tensor_2(x):
    """ Another tensor test object"""

@pytest.mark.parametrize('a', only_three_lengths)
@pytest.mark.parametrize('b', only_three_lengths)
@pytest.mark.parametrize('c', only_three_lengths)
@pytest.mark.parametrize('d', only_three_lengths)
def test_tensor_2(a, b, c, d):
    """ Test one big tensor with lots of parameters, triple self-constraint"""
    x = np.zeros((a, b, c, d))

    if a==b==c:
        tensor_2(x)
    else:
        with pytest.raises(DimensionalityError):
            tensor_2(x)
#
# Omissions
#

@dimensionality_check(y=(3,))
def omission(x,y):
    """ Only contrained on second variable"""

@pytest.mark.parametrize("a", only_three_lengths)
@pytest.mark.parametrize("b", only_three_lengths)
def test_omission(a, b):
    """ Test that we don't throw if something isn't mentioned in the decorator """
    x = np.zeros((a, ))
    y = np.zeros((b, ))

    if b == 3:
        omission(x, y)
    else:
        with pytest.raises(DimensionalityError):
            omission(x, y)
