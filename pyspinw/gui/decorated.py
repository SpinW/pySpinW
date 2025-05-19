from dataclasses import dataclass

import numpy as np

from pyspinw.site import LatticeSite

@dataclass
class InteractionFlags:
    hover: bool
    select: bool
    implied: bool

@dataclass
class DecoratedSite:
    site: LatticeSite
    position: np.ndarray
    moment: np.ndarray
    graphics_options: InteractionFlags