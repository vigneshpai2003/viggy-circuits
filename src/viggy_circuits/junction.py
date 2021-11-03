from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Set

    from .wire import Wire


class Junction:
    def __init__(self):
        self.wires: Set[Wire] = set()

    @property
    def isUsed(self) -> bool:
        """
        :return: whether the junction has been used
        """
        return len(self.wires) != 0

    @property
    def isUnused(self) -> bool:
        return not self.isUsed

    def connect(self, *wires: Wire):
        """
        not meant to be called by user, see Circuit.connect
        """
        self.wires.update(set(wires))

    def otherWires(self, *wires: Wire) -> Set[Wire]:
        """
        :returns wires other than given wires connected to this junction
        """
        return self.wires.difference(set(wires))
