from __future__ import annotations

from enum import IntEnum

from .error import CircuitError
from .switch import Switch
from .junction import Junction

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Union, Callable


class Direction(IntEnum):
    BACKWARD = -1
    NULL = 0
    FORWARD = 1


class Wire:
    def __init__(self, resistance: Union[Callable[[float], float], float],
                 battery: Union[Callable[[float], float], float, None] = None,
                 capacitance: Union[Callable[[float], float], float, None] = None,
                 inductance: Union[Callable[[float], float], float, None] = None,
                 initCharge=0.0, initCurrent=0.0):
        self.__switch = Switch()

        self.resistance = Wire.__process(resistance)
        self.battery = Wire.__process(battery)  # internal battery
        self.capacitance = Wire.__process(capacitance)
        self.inductance = Wire.__process(inductance)
        self.initCharge = initCharge
        self.initCurrent = initCurrent  # only significant if inductance is not None

        self.__junction1 = Junction()
        self.__junction2 = Junction()

    @staticmethod
    def __process(x: Union[Callable[[float], float], float, None]) -> Union[Callable[[float], float], None]:
        if callable(x) or x is None:
            return x
        else:
            return lambda t: x

    @property
    def switch(self) -> Switch:
        return self.__switch

    @property
    def junction1(self) -> Junction:
        return self.__junction1

    @property
    def junction2(self) -> Junction:
        return self.__junction2

    @property
    def junctions(self) -> Tuple[Junction, Junction]:
        return self.__junction1, self.__junction2

    @property
    def isLCR(self) -> bool:
        """
        :return: whether Wire.device has capacitor, resistor and inductor
        """
        return self.inductance is not None

    def connect(self, junction1: Junction, junction2: Junction):
        """
        not meant to be called by user, see Circuit.connect
        """
        self.__junction1 = junction1
        self.__junction2 = junction2

    @property
    def isNull(self) -> bool:
        """
        :return: true if either of the junctions has less than 2 wires or switch is open
        """
        return any([self.switch.isOpen, len(self.junction1.wires) <= 1, len(self.junction2.wires) <= 1])

    def swapJunctions(self):
        """
        swaps the internal order of junctions
        """
        self.__junction1, self.__junction2 = self.__junction2, self.__junction1

    def direction(self, junction: Junction) -> Direction:
        """
        :return: forward if junction is junction1 else backward
        """
        return Direction.FORWARD if junction is self.junction1 else Direction.BACKWARD

    def potentialDrop(self, t=0.0, q=0.0, i=0.0, di_dt=0.0):
        v = 0.0

        if self.battery is not None:
            v += self.battery(t)

        if self.resistance is not None:
            v -= i * self.resistance(t)

        if self.capacitance is not None:
            v -= q / self.capacitance(t)

        if self.inductance is not None:
            v -= di_dt * self.inductance(t)

        return v
