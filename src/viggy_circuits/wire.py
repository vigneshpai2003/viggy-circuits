from __future__ import annotations

from enum import IntEnum

from .switch import Switch
from .junction import Junction

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Union, Callable

    time_function = Callable[[float], float]
    float_like = Union[time_function, float]
    float_like_optional = Union[float_like, None]


class Direction(IntEnum):
    BACKWARD = -1
    NULL = 0
    FORWARD = 1


class Wire:
    def __init__(self, battery: float_like_optional = None, switch=None):
        self.battery = self.process(battery)
        self.__switch = switch if switch else Switch()

        self.__junction1 = Junction()
        self.__junction2 = Junction()

    @staticmethod
    def process(x: float_like_optional) -> float_like_optional:
        """
        function to process a constant and turn it into function of time
        """
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

    def potentialDrop(self, t=0.0):
        return self.battery(t) if self.battery else 0.0


class ResistorWire(Wire):
    def __init__(self, resistance: float_like,
                 battery: float_like_optional = None,
                 capacitance: float_like_optional = None,
                 initCharge=0.0, switch=None):
        super().__init__(battery, switch)

        self.resistance = Wire.process(resistance)
        self.capacitance = Wire.process(capacitance)

        self.initCharge = initCharge

        self.isLCR = False

    def potentialDrop(self, t=0.0, q=0.0, i=0.0):
        v = Wire.potentialDrop(self, t)

        if self.capacitance:
            v -= q / self.capacitance(t)

        if self.resistance:
            v -= i * self.resistance(t)

        return v


class InductorWire(ResistorWire):
    def __init__(self, inductance: float_like,
                 battery: float_like_optional = None,
                 capacitance: float_like_optional = None,
                 resistance: float_like_optional = None,
                 initCharge=0.0, initCurrent=0.0, switch=None):
        super().__init__(resistance, battery, capacitance, initCharge, switch)

        self.inductance = Wire.process(inductance)

        self.initCurrent = initCurrent

        self.isLCR = True

    def potentialDrop(self, t=0.0, q=0.0, i=0.0, di_dt=0.0):
        v = ResistorWire.potentialDrop(self, t, q, i)

        if self.inductance:
            v -= di_dt * self.inductance(t)

        return v


if TYPE_CHECKING:
    some_wire = Union[ResistorWire, InductorWire]
