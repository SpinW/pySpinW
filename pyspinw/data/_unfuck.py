""" Fix the formatting of the magnetic symmetry lookup"""

import re

# Columns should be:
# Letvin, BNS,     BNS,    OG,              OG,     INT,     INT,    MHall
# number, decimal, symbol, 3 point decimal, symbol, number, symbol, symbol

# Build a regex for that

def group(s: str):
    return r"("+s+r")"

space = r"\s+"
number = r"\d+"
point = r"\."
integer = group(number)
decimal = group(number + point + number)
double_decimal = group(number + point + number + point + number)
symbol = group(r"-?[PIFABCR].*")

regex = (integer + space +
         decimal + space +
         symbol + space +
         double_decimal + space +
         symbol + space +
         integer + space +
         symbol + space +
         symbol)

print(regex)


with open("magnetic_symmetry_conventions", 'r') as input_file:
    with open("magnetic_symmetry_conventions.csv", 'w') as output_file:
        for line in input_file:
            match = re.match(regex, line.strip())
            if match is not None:
                with_commas = ", ".join(match.groups())
                print(with_commas)
                output_file.write(with_commas)
                output_file.write("\n")