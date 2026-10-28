"""For reading and interpreting magnetic space group data."""

from dataclasses import dataclass

from importlib import resources

@dataclass
class MSGSymbols:
    """Object to hold the different kinds of symbols together.

    Parameters
    ----------
    number : int
        Sequential magnetic space group index (UNI number).
    uni : str
        Unified (UNI) magnetic space group symbol.
    bns : str
        Belov-Neronova-Smirnova (BNS) magnetic space group symbol.
    og : str
        Opechowski-Guccione (OG) magnetic space group symbol.
    """

    number: int
    uni: str
    bns: str
    og: str

def load_data():
    """Load magnetic spacegroup symbol data from package resources.

    Returns
    -------
    dict of {int: MSGSymbols}
        Mapping from magnetic space group UNI number to ``MSGSymbols`` entry.
    """
    entries = {}
    with resources.open_text("pyspinw.symmetry.data", "msg_symbols.txt") as file:
        for line in file:
            parts = [part.strip() for part in line.split("\t")]

            number = int(parts[0])

            entries[number] = MSGSymbols(number, parts[1], parts[2], parts[3])

    return entries

# To be exported
msg_symbols = load_data()
