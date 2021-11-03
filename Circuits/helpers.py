from __future__ import annotations

from math import sin
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Callable


def AC(emf_rms: float, omega: float, phi: float) -> Callable[[float], float]:
    """
    :param emf_rms:
    :param omega:
    :param phi: phase constant
    :return: function of time that returns emf
    """

    def emf(t: float) -> float:
        return 2 ** 0.5 * emf_rms * sin(omega * t + phi)

    return emf
