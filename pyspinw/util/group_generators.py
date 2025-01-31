import numpy as np
import re

from pyspinw.checks import check_sizes
from pyspinw.util.safe_expression_evaluation import evaluate_algebra

_number_regex = r"\d+(?:\.\d+)?"
_symbol_regex = r"x|y|z|\-|\+|/|\*"
_number_symbol_regex = "("+_number_regex+"|"+_symbol_regex+"|\s+)"

class Generator:

    def __init__(self,
                 rotation: np.ndarray,
                 translation: np.ndarray,
                 time_reversal: int,
                 name: str | None = None):

        self.rotation = rotation
        self.translation = translation
        self.time_reversal = time_reversal
        self._name = name

    @property
    def name(self) -> str:
        if self._name is None:
            return "<unknown>"
        else:
            return self._name

    @check_sizes(points=(-1, 6))
    def __call__(self, points: np.ndarray) -> np.ndarray:
        """ Apply this generator to a set of points and momenta"""

        new_points = points.copy()

        new_points[:, :3] = new_points[:, :3] @ self.rotation  + self.translation.reshape(1, 3)
        new_points[:, :3] %= 1 # To unit cell

        new_points[:, 3:] *= self.time_reversal

        return new_points

    def and_then(self, other: "Generator") -> "Generator":
        """ Composition of generators """

        return Generator(
            rotation = self.rotation @ other.rotation,
            translation = (self.translation.reshape(1, 3) @ other.rotation + other.translation) % 1,
            time_reversal = self.time_reversal * other.time_reversal,
            name=self.name + "->" + other.name)

    def __lt__(self, other: "Generator") -> bool:
        raise NotImplementedError()

def _convert_token_to_number(token: str, to_zero: str, to_one: str) -> str:
    """ Helper function, converts a token to "0" or "1" based on whether it is in each of the given strings

    e.g. _convert_token_to_number("y", "xy", "z") yields "0"
         _convert_token_to_number("z", "xy", "z") yields "1"

    :param token: the token to be potentially converted
    :param to_zero: tokens in this string will yield "0", this has priority
    :param to_one: tokens in this string will yield "1"
    :returns: "0", "1" or the input token
    """

    if token in to_zero:
        return "0"
    elif token in to_one:
        return "1"
    else:
        return token

def _evaluate_generator_with_subsitution(tokens: list[str], to_zero: str, to_one: str):
    """ Helper function: Evaluate a tokenised (see parse_space_group_generator) list


    :param tokens: list of tokens representing the generator
    :param to_zero: tokens in this string will become zeros
    :param to_one: tokens in this string will become ones
    :returns: the value of the function evaluated with the specified substitutions
    """

    with_numbers = "".join([_convert_token_to_number(token, to_zero, to_one) for token in tokens])
    return evaluate_algebra(with_numbers)

def parse_space_group_generator(
        generator_string: str,
        time_reversed: bool | None = None) -> tuple[np.ndarray, np.ndarray, float | None]:

    """ Parse a space group generator string, e.g. '-x,y,-z+1/2' (three components)
    or 'x+1/2,y+1/2,z,-1' (four components, magnetic)

    :returns: 'rotation' matrix, translation, and time reversal
    """

    components = generator_string.split(",")

    if time_reversed is None:
        if len(components) == 3:
            time_reversal = 1.0
        elif len(components) == 4:
            time_reversal = float(components[3])
        else:
            raise ValueError("Expected three or four comma separated values")

    else:
        if len(components) != 3:
            raise ValueError("Expected exactly three comma separated values for case with time reversal specified")

        time_reversal = -1.0 if time_reversed else 1.0

    # parse the linear equations specifying the magnetic space group
    # general strategy is to just use pythons abstract syntax tree to evaluate the linear
    # expression at different x,y,z values to deduce the constants in a.(x,y,z) + b for each one,
    # then build the appropriate matrix

    quadratic = []
    linear = []

    # Tokenise to check its sanitary
    for component in components:
        tokens = re.findall(_number_symbol_regex, component)
        sanitised = "".join(tokens)

        if component != sanitised:
            raise ValueError(f"Invalid generator string ('{component}' does not match sanitised '{sanitised}')")

        # Find constant by substituting 0 for x,y and z
        b = _evaluate_generator_with_subsitution(tokens, "xyz", "")

        # Same for each of x,y and z, but subtracting the constant
        a_x = _evaluate_generator_with_subsitution(tokens, "yz", "x") - b
        a_y = _evaluate_generator_with_subsitution(tokens, "xz", "y") - b
        a_z = _evaluate_generator_with_subsitution(tokens, "xy", "z") - b

        a = [a_x, a_y, a_z]

        quadratic.append(a)
        linear.append(b)

    quadratic = np.array(quadratic)
    linear = np.array(linear)

    return quadratic, linear, time_reversal

def parse_one_line_generators(generator_string: str):
    """ Expects the generators to be a single line of tuples in terms of x,y,z, separated by commans

    e.g. (-x,y,-z+1/2);(x,-y,z+1/2);(x+1/2,y+1/2,z) """

    individual_strings = generator_string.split(";")
    output = []

    for string in individual_strings:

        if string.endswith("'"):
            time_reversed = True
            string = string[:-1]
        else:
            time_reversed = False

        if string.startswith("(") and string.endswith(")"):
            string = string[1:-1]
        else:
            raise ValueError("Expected brackets around generator strings")


        output.append(parse_space_group_generator(string, time_reversed=time_reversed))

    return output

def spglib_generators_to_list(generators: dict) -> list[tuple[np.ndarray, np.ndarray, float]]:
    rotations = generators["rotations"]
    translations = generators["translations"]
    time_reversals = generators["time_reversals"]

    return [(rotations[i,:,:], translations[i,:], -1.0 if time_reversals[i] > 0.5 else 1.0)
            for i in range(len(time_reversals))]

if __name__ == "__main__":
    from pyspinw.util.magnetic_symmetry import name_converter

    test_string = name_converter.litvin[97].generators

    generators = parse_one_line_generators(test_string)

    print(generators)