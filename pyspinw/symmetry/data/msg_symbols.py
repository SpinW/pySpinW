from dataclasses import dataclass

from importlib import resources

@dataclass
class MSGSymbols:
    """ Object to hold the different kinds of symbols together"""
    number: int
    uni: str
    bns: str
    og: str

def load_data():
    """ Main method for loading magnetic spacegroup data"""
    entries = {}
    with resources.open_text("pyspinw.symmetry.data", "msg_symbols.txt") as file:
        for line in file:
            parts = [part.strip() for part in line.split("\t")]

            number = int(parts[0])

            entries[number] = MSGSymbols(number, parts[1], parts[2], parts[3])

    return entries

# To be exported
msg_symbols = load_data()