""" Classes containing extra data used for rendering """

from dataclasses import dataclass

import numpy as np

from pyspinw.site import LatticeSite

@dataclass
class InteractionFlags:
    """ State of the user interaction with this site"""

    hover: bool
    selected: bool
    implied: bool
    implied_hover: bool

@dataclass
class DecoratedSite:
    """ Site with extra information for drawing """

    site: LatticeSite
    position: np.ndarray
    moment: np.ndarray
    flags: InteractionFlags
