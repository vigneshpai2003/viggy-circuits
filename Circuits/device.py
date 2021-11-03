from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Union, Callable


class Device:
    def __init__(self, resistance: Union[Callable[[float], float], float],
                 battery: Union[Callable[[float], float], float, None] = None,
                 capacitance: Union[Callable[[float], float], float, None] = None,
                 inductance: Union[Callable[[float], float], float, None] = None,
                 initCharge=0.0, initCurrent=0.0):
        """
        :param resistance: a float or function of time that returns float
        :param battery: optional
        :param capacitance: optional
        :param inductance: optional
        :param initCharge: optional, only need to change if capacitance is given
        :param initCurrent: optional, only need to change if inductance is given
        """
        self.resistance = Device.__process(resistance)
        self.battery = Device.__process(battery)  # internal battery
        self.capacitance = Device.__process(capacitance)
        self.inductance = Device.__process(inductance)
        self.initCharge = initCharge
        self.initCurrent = initCurrent  # only significant if inductance is not None

    @staticmethod
    def __process(x: Union[Callable[[float], float], float, None]) -> Union[Callable[[float], float], None]:
        if callable(x) or x is None:
            return x
        else:
            return lambda t: x

    def potentialDrop(self, t: float, q=0.0, dq_dt=0.0, d2q_dt2=0.0):

        v = 0.0

        if self.battery is not None:
            v += self.battery(t)

        if self.resistance is not None:
            v -= dq_dt * self.resistance(t)

        if self.capacitance is not None:
            v -= q / self.capacitance(t)

        if self.inductance is not None:
            v -= d2q_dt2 * self.inductance(t)

        return v
