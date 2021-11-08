from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Set

    from .wire import some_wire


class Junction:
    def __init__(self):
        self.wires: Set[some_wire] = set()

    @property
    def isUsed(self) -> bool:
        """
        :return: whether the junction has been used
        """
        return len(self.wires) != 0

    def connect(self, wire: some_wire):
        """
        not meant to be called by user, see Circuit.connect
        """
        self.wires.add(wire)
