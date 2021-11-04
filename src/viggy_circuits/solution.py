from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from typing import Dict, List, Iterable
    from .wire import Wire


class WireSolution:
    def __init__(self, wire: Wire):
        self.__wire = wire
        self.t = []
        self.q = []
        self.i = []
        self.di_dt = [] if wire.isLCR else None

    def freeze(self):
        self.t = np.array(self.t)
        self.q = np.array(self.q)
        self.i = np.array(self.i)
        if self.__wire.isLCR:
            self.di_dt = np.array(self.di_dt)

        # potential difference across wire ends
        self.v = np.frompyfunc(self.__wire.device.potentialDrop, 4, 1)(self.t, self.q, self.i, self.di_dt)

        # power dissipated across resistor
        self.powerR = np.square(self.i) * np.vectorize(self.__wire.device.resistance)(self.t)

        # energy stored in capacitor
        self.energyC = 0.5 * np.square(self.q) / np.vectorize(self.__wire.device.capacitance)(self.t)

        # energy stored in inductor
        if self.__wire.isLCR:
            self.energyL = 0.5 * np.square(self.i) * np.vectorize(self.__wire.device.inductance)(self.t)

    def __iter__(self):  # allows unpacking of this object
        """
        :return: a tuple of (t, q, i, di_dt) lists, di_dt is absent for non inductor wire
        """
        return ([self.t, self.q, self.i] + [self.di_dt] if self.__wire.isLCR else []).__iter__()

    @property
    def lastQ(self):
        return self.q[-1] if self.q else self.__wire.device.initCharge

    @property
    def lastI(self):
        return self.i[-1] if self.i else self.__wire.device.initCurrent

    def update(self, t, q, i, di_dt=None):
        self.t.append(t)
        self.q.append(q)
        self.i.append(i)
        if di_dt is not None:
            self.di_dt.append(di_dt)

    def updateVoid(self, t):
        self.t.append(t)
        self.q.append(self.lastQ)
        self.i.append(0.0)
        if self.__wire.isLCR:
            self.di_dt.append(0.0)


class Solution:
    def __init__(self, wires: Iterable[Wire]):
        self.__solutionOf: Dict[Wire, WireSolution] = dict([(wire, WireSolution(wire)) for wire in wires])

    def __getitem__(self, item: Wire):
        return self.__solutionOf[item]

    def freeze(self):
        for sol in self.__solutionOf.values():
            sol.freeze()

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
