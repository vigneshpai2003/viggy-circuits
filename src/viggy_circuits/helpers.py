from __future__ import annotations

import numpy as np
from math import tau
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Callable


def generatorAC(emf: float, f: float, phi: float) -> Callable[[float], float]:
    """
    :param emf: root mean squared value of emf
    :param f: frequency
    :param phi: phase constant
    :return: function of time that returns emf
    """

    def _AC(t: float) -> float:
        return np.sqrt(2) * emf * np.sin(tau * f * t + phi)

    return _AC
