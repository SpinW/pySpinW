from dataclasses import dataclass

import numpy as np

from pyspinw.site import LatticeSite

@dataclass
class InteractionFlags:
    hover: bool
    selected: bool
    implied: bool
    implied_hover: bool

@dataclass
class DecoratedSite:
    site: LatticeSite
    position: np.ndarray
    moment: np.ndarray
    flags: InteractionFlags
