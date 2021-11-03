from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Set

    from .wire import Wire


class Junction:
    def __init__(self):
        self.wires: Set[Wire] = set()

    @property
    def is_used(self) -> bool:
        return len(self.wires) != 0.0

    @property
    def is_unused(self) -> bool:
        return not self.is_used

    def connect(self, *wires: Wire):
        """
        not meant to be called directly, see Circuit.connect
        """
        self.wires.update(set(wires))

    def other_wires(self, *wires: Wire) -> Set[Wire]:
        """
        :returns wires other than given wires connected to this junction
        """
        return self.wires.difference(set(wires))
