from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Union, Callable


class Device:
    def __init__(self, resistance: Union[Callable[[float], float], float],
                 battery: Union[Callable[[float], float], float, None] = None,
                 capacitance: Union[Callable[[float], float], float, None] = None,
                 inductance: Union[Callable[[float], float], float, None] = None,
                 init_charge=0.0, init_current=0.0):
        """
        :param resistance: a float or function of time that returns float
        :param battery: optional
        :param capacitance: optional
        :param inductance: optional
        :param init_charge: optional, only need to change if capacitance is given
        :param init_current: optional, only need to change if inductance is given
        """
        self.resistance = Device.process(resistance)
        self.battery = Device.process(battery)  # internal battery
        self.capacitance = Device.process(capacitance)
        self.inductance = Device.process(inductance)
        self.init_charge = init_charge
        self.init_current = init_current  # only significant if inductance is not None

    @classmethod
    def process(cls, x: Union[Callable[[float], float], float, None]) -> Union[Callable[[float], float], None]:
        if callable(x) or x is None:
            return x
        else:
            return lambda t: x

    def potential_drop(self, t: float, q: float = 0.0, dq_dt: float = 0.0, d2q_dt2: float = 0.0):

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
