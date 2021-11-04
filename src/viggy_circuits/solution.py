from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Dict, List, Iterable
    from .wire import Wire


class WireSolution:
    def __init__(self, wire: Wire):
        self.__wire = wire
        self.__t: List[float] = []
        self.__q: List[float] = []
        self.__i: List[float] = []
        self.__di_dt: List[float] = [] if wire.isLCR else None

    def __iter__(self):  # allows unpacking of this object
        """
        :return: a tuple of (t, q, i, di_dt) lists, di_dt is absent for non inductor wire
        """
        return ([self.__t, self.__q, self.__i] + [self.__di_dt] if self.__wire.isLCR else []).__iter__()

    @property
    def lastQ(self):
        return self.__q[-1] if self.__q else self.__wire.device.initCharge

    @property
    def lastI(self):
        return self.__i[-1] if self.__i else self.__wire.device.initCurrent

    def update(self, t, q, i, di_dt=None):
        self.__t.append(t)
        self.__q.append(q)
        self.__i.append(i)
        if di_dt is not None:
            self.__di_dt.append(di_dt)

    def updateVoid(self, t):
        self.__t.append(t)
        self.__q.append(self.lastQ)
        self.__i.append(0.0)
        if self.__wire.isLCR:
            self.__di_dt.append(0.0)


class Solution:
    def __init__(self, wires: Iterable[Wire]):
        self.__solutionOf: Dict[Wire, WireSolution] = dict([(wire, WireSolution(wire)) for wire in wires])

    def __getitem__(self, item: Wire):
        return self.__solutionOf[item]

    def lastKnownValues(self, wires: List[Wire]):
        last = []
        for wire in wires:
            last.append(self[wire].lastQ)
            if wire.isLCR:
                last.append(self[wire].lastI)

        return last

    @staticmethod
    def mapWireToIndex(wires: List[Wire]):
        # map each wire to index of corresponding value in lastKnownValues
        indexMap = {}

        i = 0
        for wire in wires:
            indexMap[wire] = i
            i += 1
            if wire.isLCR:  # LCR wires have two corresponding values
                i += 1

        return indexMap
