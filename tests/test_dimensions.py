""" Tests for the check_dimensions decorator

Soooo many decorators, it must be Christmas.

It's November.

So yes.
"""
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
    input = np.zeros((length,))

    # Check vector with right name
    if not success:
        with pytest.raises(DimensionalityError):
            fixed_length_vector(input)

        # Should work if keyword arguments are supplied instead
        with pytest.raises(DimensionalityError):
            fixed_length_vector(x=input)

    else:

        # TODO: check this is how to test a function just runs without error
        fixed_length_vector(input)
        fixed_length_vector(x=input)

    # The following this should throw
    #  * wrong variable name specified
    #  * wrong dimensionality

    with pytest.raises(Exception):
        fixed_length_vector(y=input)

    with pytest.raises(DimensionalityError):
        fixed_length_vector(np.zeros((4,5)))


@pytest.mark.parametrize("function", [string_arbitrary_length_vector,
                                      int_arbitrary_length_vector])
@pytest.mark.parametrize("length", vector_check_lengths)
def test_vector_parameter(function, length):
    input = np.zeros((length,))

    # Should just not fail
    function(input)
    function(x=input)

    # The following this should throw
    #  * wrong variable name specified
    #  * wrong dimensionality

    with pytest.raises(Exception):
        function(y=input)

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
    mat = np.zeros((a, b))
    function(mat)

@dimensionality_check(x=(1,7))
def one_by_seven_matrix(x):
    pass

@pytest.mark.parametrize('a', vector_check_lengths)
@pytest.mark.parametrize('b', vector_check_lengths)
def test_fixed_size_matrix(a, b):
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

    should_pass = xa == zb and xb == yb and ya == za

    x = np.zeros((xa, xb))
    y = np.zeros((ya, yb))
    z = np.zeros((za, zb))

    if should_pass:
        triple_constraint(x,y,z)
    else:
        with pytest.raises(DimensionalityError):
            triple_constraint(x,y,z)


#
# Omissions
#