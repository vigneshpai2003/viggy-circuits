from __future__ import annotations

from .error import CircuitError
from .switch import Switch
from .junction import Junction

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple

    from .device import Device


class Wire:
    def __init__(self, device: Device, junction1: Junction = None, junction2: Junction = None):

        self.__switch = Switch()

        self.__device = device

        if junction1 is None:
            self.__junction1 = Junction()
        else:
            self.__junction1 = junction1

        if junction2 is None:
            self.__junction2 = Junction()
        else:
            self.__junction2 = junction1

    @property
    def device(self) -> Device:
        return self.__device

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

    def is_LCR(self) -> bool:
        """
        whether wire.device has capacitor, resistor and inductor
        """
        return self.__device.inductance is not None

    def connect(self, junction: Junction):
        """
        not meant to be called directly, see Circuit.connect
        """
        if self.connection_complete:
            raise CircuitError(f"Cannot connect 3rd wire to {self}")

        if len(self.__junction1.wires) == 0:
            self.__junction1 = junction
        elif len(self.__junction2.wires) == 0:
            self.__junction2 = junction

    @property
    def connection_complete(self) -> bool:
        return self.__junction1.is_used and self.__junction2.is_used

    @property
    def is_null(self) -> bool:
        return self.__switch.is_open or not self.connection_complete \
               or len(self.__junction1.wires) == 1 or len(self.__junction2.wires) == 1

    def swap_junctions(self):
        self.__junction1, self.__junction2 = self.__junction2, self.__junction1

    def other_junction(self, junction: Junction):
        if junction is self.__junction1:
            return self.__junction2
        elif junction is self.__junction2:
            return self.__junction1
        else:
            raise CircuitError(f"given Junction: {junction} is not connected to wire")
