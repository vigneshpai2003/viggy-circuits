from __future__ import annotations

import copy
import random

import numpy as np
from scipy import linalg, integrate

from .error import CircuitError
from .wireCollection import WireCollection, Direction

from typing import TYPE_CHECKING, FrozenSet

if TYPE_CHECKING:
    from typing import Set, List, Dict, Tuple, Callable, Any, Any

    from .wire import Wire
    from .junction import Junction


class Circuit:
    def __init__(self):
        # all the wires and junctions in the circuit
        self.__wires: Set[Wire] = set()
        self.__junctions: Set[Junction] = set()

        # all the wires and junctions that are required to simulate (changes every iteration)
        self.__effectiveWires: Set[Wire] = set()
        self.__effectiveJunctions: Set[Junction] = set()

        # all the wires and junctions that do not contribite to circuit (changes every iteration)
        self.__voidWires: Set[Wire] = set()
        self.__voidJunctions: Set[Junction] = set()

        self.__events = {}

    def connect(self, junction: Junction, *wires: Wire):
        """
        function that connects junction to wire and each wire to junction
        """
        self.__junctions.add(junction)
        self.__wires.update(wires)

        junction.connect(*wires)
        for wire in wires:
            wire.connect(junction)

    def __simplifyCircuit(self):
        """
        sorts the wires and junctions as effective or void
        """
        # reset the last simplification
        self.__voidWires = set()
        self.__voidJunctions = set()

        # all null wires are automatically void
        for wire in self.__wires:
            if wire.isNull:
                self.__voidWires.add(wire)

        # make it -1 so it enters while loop
        voidJunctionCount = -1

        while voidJunctionCount != 0:  # if no void junctions are found
            voidJunctionCount = 0

            # find void junctions
            for junction in self.__junctions - self.__voidJunctions:
                # (junction.wires - self.__voidWires) is the effective wires attached to the junction
                if len(junction.wires - self.__voidWires) in (0, 1):
                    voidJunctionCount += 1
                    self.__voidJunctions.add(junction)

                    # if junction is void, all its wires also become void
                    self.__voidWires.update(junction.wires)

        self.__effectiveWires = self.__wires - self.__voidWires
        self.__effectiveJunctions = self.__junctions - self.__voidJunctions

    def __getLoopsWith(self, wire: Wire) -> List[WireCollection]:
        """
        this function should be called after simplifying circuit
        :param wire: should be in self.__effectiveWires
        :return: all the closed loops with the wire
        """
        loops = []
        endpoint = wire.junction1

        def __loops(currentJunction: Junction, currentWire: Wire,
                    loopPath: WireCollection, parsedJunctions: Set[Junction]):
            """
            :param currentJunction: the junction up to which the path was calculated
            :param currentWire: the wire up to which the path was calculated
            :param loopPath: the path of the loop so far
            :param parsedJunctions: the junctions that are already a part of the loop
            :return:
            """
            for newWire in currentJunction.otherWires(currentWire) - self.__voidWires:
                """
                for every wire in the junction, an attempt is made to add it to loop_path
                if failed, the entire loop path corresponding up to that point is discarded
                """
                # direction is assumed to be from junction1 to junction2 always
                if currentJunction is newWire.junction1:
                    sign = Direction.FORWARD
                elif currentJunction is newWire.junction2:
                    sign = Direction.BACKWARD
                else:
                    raise CircuitError("sign of wire could not be deduced")

                # newLoopPath must be copied and not referenced
                newLoopPath = copy.copy(loopPath)
                newLoopPath.append(newWire, sign)

                newJunction = newWire.otherJunction(currentJunction)

                if newJunction is endpoint:  # if closed loop is found
                    loops.append(newLoopPath)
                    continue
                elif newJunction in parsedJunctions:  # if junction was already encountered before
                    continue

                __loops(newJunction, newWire, newLoopPath, parsedJunctions.union({newJunction}))

        __loops(currentJunction=wire.junction2, currentWire=wire,
                loopPath=WireCollection((wire, Direction.FORWARD)), parsedJunctions={wire.junction2})

        return loops

    def __firstLawWires(self) -> Set[WireCollection]:
        """
        :return: the wires to be used in implementation of Kirchhoff's First Law
        """
        equations: Set[WireCollection] = set()

        # possibleEquations is updated to include all possible combinations of already formed equations
        possibleEquations: Set[FrozenSet[Wire]] = set()

        for junction in self.__effectiveJunctions:
            # remove wires that are void
            junctionWires = junction.wires - self.__voidWires

            if all([wire.isLCR for wire in junctionWires]):  # if all wires have inductor
                effectiveJunctionWires = junctionWires.copy()
            else:
                # if some or no wires have inductor
                # current in all wires with inductors are known in __derivatives
                # only those without inductor should be included
                effectiveJunctionWires = set()
                for wire in junctionWires:
                    if not wire.isLCR:
                        effectiveJunctionWires.add(wire)

            # check if these wire can be derived from existing wire equations
            if effectiveJunctionWires in possibleEquations:
                continue

            # update possible combos
            possibleEquations.add(frozenset(effectiveJunctionWires))  # cant create set of sets
            possibleEquations.update(set(equationWires.symmetric_difference(effectiveJunctionWires)
                                         for equationWires in possibleEquations))

            """
            if 2 equation in the currents of a certain set of wire exists,
            another equation can be created by subtracting them,
            this new equation will not contain variables common to both equations,
            it will contain the symmetric difference of the variables
            Note: the orientation of wire is not being considered, should it be?
            """

            loop = WireCollection()
            for wire in junctionWires:
                if junction is wire.junction1:
                    loop.append(wire, Direction.FORWARD)
                else:
                    loop.append(wire, Direction.BACKWARD)
            equations.add(loop)

        return equations

    def __secondLawWires(self, limit=None) -> Set[WireCollection]:
        """
        :param limit: the number of equations needed
        return the wires to be used in implementation of Kirchhoff's Second Law;
        """
        equations = set()
        parsedWires = set()

        for wire in self.__effectiveWires:
            for loop in self.__getLoopsWith(wire):
                if set(loop.wires).issubset(parsedWires):  # if equation in these variablles is already formed
                    continue

                parsedWires.update(set(loop.wires))
                equations.add(loop)

                if len(equations) == limit:  # if limit is reached no need to do further iterations
                    return equations

        return equations

    @staticmethod
    def __derivatives(t: float, x: np.ndarray, wires: List[Wire], indices: Dict[Wire, int],
                      firstLawEquations: Set[WireCollection],
                      secondLawEquations: Set[WireCollection]) -> np.ndarray:
        """
        :param t: time
        :param x: a list of charge and current values
        :param wires: the list of wires corresponding to each value in x
        :param indices: the index of the charge / current value in x for each wire

        x[indices[wire]] gives charge
        x[indices[wire] + 1] gives current (only applicable for LCR)

        :return: the corresponding derivatives of x
        """
        degree = len(wires)  # number of variables to solve for

        # R i = V
        R = np.zeros((degree, degree))
        V = np.zeros(degree)

        # fill R, V using first law
        i = 0
        for equation in firstLawEquations:
            allLCR = all([wire.isLCR for wire in equation.wires])
            for wire in equation.wires:
                sign = equation.getSign(wire)
                if wire.isLCR and not allLCR:
                    V[i] -= sign * x[indices[wire] + 1]
                else:
                    R[i, wires.index(wire)] = sign

            i += 1

        # fill R, V using second law
        for equation in secondLawEquations:
            for wire in equation.wires:
                sign = equation.getSign(wire)

                if wire.isLCR:
                    R[i, wires.index(wire)] = sign * wire.device.inductance(t)

                    # calculate potential drop due to resistor
                    V[i] -= sign * wire.device.resistance(t) * x[indices[wire] + 1]
                else:
                    R[i, wires.index(wire)] = sign * wire.device.resistance(t)

                # potential drop due to battery
                if wire.device.battery is not None:
                    V[i] += sign * wire.device.battery(t)

                # potential drop due to capacitor
                if wire.device.capacitance is not None:
                    V[i] -= sign * x[indices[wire]] / wire.device.capacitance(t)

            i += 1

        derivatives = linalg.solve(R, V)

        dx_dt = np.zeros_like(x)

        for wire in wires:
            if wire.isLCR:
                dx_dt[indices[wire]] = x[indices[wire] + 1]  # i is already in input
                dx_dt[indices[wire] + 1] = derivatives[wires.index(wire)]  # di_dt
            else:
                dx_dt[indices[wire]] = derivatives[wires.index(wire)]  # i

        return dx_dt

    @staticmethod
    def __indexMap(wires: List[Wire]) -> Dict[Wire, int]:
        # map each wire to index of corresponding value in initConditions
        indexMap = {}

        i = 0
        for wire in wires:
            indexMap[wire] = i
            i += 1
            if wire.isLCR:  # LCR wires have two corresponding values
                i += 1

        return indexMap

    def solve(self, end: float, dt: float) -> Dict[Wire, Tuple[List[float], List[List[float]]]]:
        # empty initialization of solution dictionary
        solution = dict()
        for wire in self.__wires:
            solution[wire] = ([], [])

        # calculate effective wires and junctions
        self.__simplifyCircuit()
        degree = len(self.__effectiveWires)
        wires = random.sample(list(self.__effectiveWires), degree)

        # set initial conditions of simulation
        initConditions = []
        for wire in wires:
            initConditions.append(wire.device.initCharge)
            if wire.isLCR:
                initConditions.append(wire.device.initCurrent)

        for (start, stop), events in self.__getIntervals(end):
            indexMap = Circuit.__indexMap(wires)

            # calculate firstLaw and secondLaw
            firstLaw = self.__firstLawWires()
            secondLaw = self.__secondLawWires(limit=degree - len(firstLaw))

            # method: Nonstiff = ('RK45', 'RK23', 'DOP853'); Stiff = ('Radau', 'BDF'); Universal = 'LSODA'
            integrated_sol = integrate.solve_ivp(fun=Circuit.__derivatives, y0=initConditions,
                                                 t_span=(start, stop), t_eval=np.arange(start, stop, dt),
                                                 args=(wires, indexMap, firstLaw, secondLaw))

            # update solution
            for t_index in range(len(integrated_sol.t)):
                time = integrated_sol.t[t_index]

                x = integrated_sol.y[:, t_index]

                derivatives = self.__derivatives(time, x, wires, indexMap, firstLaw, secondLaw)

                # add all values to solution
                for wire in wires:
                    time_list, q_list = solution[wire]
                    time_list.append(time)
                    q = integrated_sol.y[indexMap[wire], t_index]

                    if wire.isLCR:
                        i = integrated_sol.y[indexMap[wire] + 1, t_index]
                        di_dt = derivatives[indexMap[wire] + 1]
                        q_list.append((q, i, di_dt))
                    else:
                        i = derivatives[indexMap[wire]]
                        q_list.append((q, i))

                # current and its derivative is 0 for void_wires, charge is init_charge or last known value
                for wire in self.__voidWires:
                    time_list, q_list = solution[wire]
                    time_list.append(time)

                    if len(q_list) == 0:
                        charge = wire.device.initCharge
                    else:
                        charge = q_list[-1][0]

                    if wire.isLCR:
                        q_list.append((charge, 0.0, 0.0))
                    else:
                        q_list.append((charge, 0.0))

            # call event functions
            for func, args, kwargs in events:
                func(*args, **kwargs)

            if not stop == end:  # if this is not the last iteration
                # simplify circuit again
                self.__simplifyCircuit()
                degree = len(self.__effectiveWires)
                wires = random.sample(list(self.__effectiveWires), degree)

                # set initConditions to the last known values
                initConditions = []
                for wire in wires:
                    time_list, q_list = solution[wire]
                    initConditions.extend(q_list[-1][:-1])

        return solution

    def addEvent(self, time: float, event: Callable, args=None, kwargs=None):
        """
        when solving, the event is called using args and kwargs at specified time;
        multiple events are allowed at same time, but add_event must be called multiple times;
        """
        if args is None:
            args = tuple()
        if kwargs is None:
            kwargs = dict()

        if time in self.__events.keys():
            self.__events[time].append((event, args, kwargs))
        else:
            self.__events[time] = [(event, args, kwargs)]

    def __getIntervals(self, end: float) -> List[Tuple[Tuple[float, float], List[Tuple[Callable, Any, Any]]]]:
        """
        event = {t0 = [events, ...], t1 = [events, ...], ...};
        splits time interval into smaller intervals separated by event execution;
        """
        if len(self.__events) == 0:
            return [((0.0, end), [])]

        times = sorted([i for i in list(self.__events.keys()) if i <= end])
        intervals = []
        old = 0.0
        for i in range(len(times)):
            intervals.append(((old, times[i]), self.__events[times[i]]))
            old = times[i]

        intervals.append(((times[-1], end), []))

        return intervals
