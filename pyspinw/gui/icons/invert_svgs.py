""" Script to make darkmode version of icons """

import os
import re

# Matches hex colours that are repeated
hex_color_pattern = re.compile(r'#([0-9a-fA-F]{2})\1\1')

def invert_gray_match(m):
    """ Invert an individual match"""
    value = int(m.group(1), 16)
    inverted = 255 - value
    hex_byte = f'{inverted:02x}'
    return f'#{hex_byte * 3}'


directory = "svg"

for filename in os.listdir(directory):


    if filename.endswith(".svg") and not filename.endswith("-dark.svg"):
        print("Translating", filename)
        with open(os.path.join(directory, filename), 'r') as fin:
            data = fin.read()

            with open(os.path.join(directory, filename[:-4] +  "-dark.svg"), 'w') as fout:
                # we want to match anything of the form
                #  "#PQPQPQ" and replace it with "#RSRSRS" where
                # RS + PQ = 0xFF

                fout.write(hex_color_pattern.sub(invert_gray_match, data))

    else:
        print("Skipping", filename)
