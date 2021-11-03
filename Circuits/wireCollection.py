from __future__ import annotations

from enum import IntEnum

from .error import CircuitError

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, List

    from .wire import Wire


class Direction(IntEnum):
    BACKWARD = -1
    NULL = 0
    FORWARD = 1


class WireCollection:
    def __init__(self, *args: Tuple[Wire, Direction]):
        self.__wires: List[Wire] = []
        self.__directions: List[Direction] = []
        for arg in args:
            self.append(*arg)

    @property
    def wires(self) -> List[Wire]:
        return self.__wires

    @property
    def directions(self) -> List[Direction]:
        return self.__directions

    def __iter__(self):
        self.__n = -1
        return self

    def __next__(self):
        if self.__n < len(self) - 1:
            self.__n += 1
            return self.__wires[self.__n], self.__directions[self.__n]
        else:
            raise StopIteration

    def __len__(self):
        return len(self.__wires)

    def __copy__(self) -> WireCollection:
        new_obj = WireCollection()
        new_obj.__wires = [wire for wire in self.__wires]
        new_obj.__directions = [direction for direction in self.__directions]
        return new_obj

    def append(self, wire: Wire, direction: Direction):
        if wire in self.__wires:
            raise CircuitError(f"Wire: {wire} already exists in SignedWires: {self}")
        self.__wires.append(wire)
        self.__directions.append(direction)

    def getDirection(self, wire: Wire) -> Direction:
        return self.__directions[self.__wires.index(wire)]
