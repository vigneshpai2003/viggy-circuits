from __future__ import annotations

import copy
import random

from numpy import arange
from scipy import linalg, integrate

from .error import CircuitError
from .wireCollection import WireCollection, Direction

from typing import TYPE_CHECKING, FrozenSet

if TYPE_CHECKING:
    from typing import Set, List, Dict, Union, Tuple, Callable, Any, Any

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

    def __derivatives(self, t: float, q: Dict[Wire, Union[Tuple[float], Tuple[float, float]]],
                      first_law_equations: Set[WireCollection],
                      second_law_equations: Set[WireCollection]) -> Dict[Wire, float]:
        """
        t is time and q is dictionary corresponding each wire to their respective charges

        q[wire] is of the form [q,] for wires without inductor
        q[wire] is of the form [q, i] for wires with inductor

        returns dq_dt for wires without inductor and d2q_dt2 for wires with inductor in the form of a dictionary
        """
        degree = len(self.__effectiveWires)  # number of variables to solve for
        wires = random.sample(tuple(self.__effectiveWires), degree)  # this will indicate order of wires in matrix

        # R i = V
        R = []
        V = []

        # fill R, V using first law
        for equation in first_law_equations:
            R_eq = [0.0] * degree
            V_eq = 0.0

            if all([wire.isLCR for wire in equation.wires]):  # all are LCR
                for wire in equation.wires:
                    R_eq[wires.index(wire)] = float(equation.getDirection(wire))
            else:
                for wire in equation.wires:
                    if wire.isLCR:  # if LCR then dq_dt is known
                        V_eq -= float(equation.getDirection(wire)) * q[wire][1]
                    else:
                        R_eq[wires.index(wire)] = float(equation.getDirection(wire))

            R.append(R_eq)
            V.append(V_eq)

        # fill R, V using second law
        for equation in second_law_equations:
            R_eq = [0.0] * degree
            V_eq = 0.0

            for wire in equation.wires:
                sign = float(equation.getDirection(wire))

                if wire.isLCR:  # if LCR
                    R_eq[wires.index(wire)] = sign * wire.device.inductance(t)

                    # calculate potential drop due to resistor
                    V_eq -= sign * wire.device.resistance(t) * q[wire][1]
                else:
                    R_eq[wires.index(wire)] = sign * wire.device.resistance(t)

                # potential drop due to battery and capacitor
                if wire.device.battery is not None:
                    V_eq += sign * wire.device.battery(t)

                if wire.device.capacitance is not None:
                    V_eq -= sign * q[wire][0] / wire.device.capacitance(t)

            R.append(R_eq)
            V.append(V_eq)

        solution = linalg.solve(R, V)

        derivatives = dict()

        for i in range(degree):
            derivatives[wires[i]] = solution[i]

        return derivatives

    def solve(self, end: float, dt: float) -> Dict[Wire, Tuple[List[float], List[List[float]]]]:
        # empty initialization of solution dictionary
        solution = dict()
        for wire in self.__wires:
            solution[wire] = ([], [])

        # calculate effective wires and junctions
        self.__simplifyCircuit()
        degree = len(self.__effectiveWires)
        wires = random.sample(list(self.__effectiveWires), degree)

        # set initial conditions
        init_conditions = []
        for wire in wires:
            if wire.device.inductance is not None:
                # if LCR
                init_conditions.extend([wire.device.initCharge, wire.device.initCurrent])
            else:
                # if RC
                init_conditions.append(wire.device.initCharge)

        for event in self.__get_intervals(end):
            # calculate first_law and second_law
            first_law = self.__firstLawWires()
            second_law = self.__secondLawWires(limit=degree - len(first_law))

            def derivative_wrapper(_t: float, _q: List[float]):
                # convert x as a list to a dict to pass to self.__derivatives
                q_dict = {}
                _i = 0
                for _wire in wires:
                    if _wire.device.inductance is not None:
                        # if LCR
                        q_dict[_wire] = (_q[_i], _q[_i + 1])
                        _i += 2
                    else:
                        # if RC
                        q_dict[_wire] = (_q[_i],)
                        _i += 1

                derivatives_dict = self.__derivatives(_t, q_dict, first_law, second_law)

                # convert dict obtained back to list with derivatives of q
                _derivatives = []
                for _wire in wires:
                    if _wire.device.inductance is not None:
                        # if LCR
                        _derivatives.extend([q_dict[_wire][1], derivatives_dict[_wire]])
                    else:
                        # if RC
                        _derivatives.append(derivatives_dict[_wire])
                return _derivatives

            # method = Nonstiff = ('RK45', 'RK23', 'DOP853'); Stiff = ('Radau', 'BDF'); Universal = 'LSODA'
            integrated_sol = integrate.solve_ivp(fun=derivative_wrapper, y0=init_conditions,
                                                 t_span=event[0], t_eval=arange(*event[0], dt))

            # update solution
            for t_index in range(len(integrated_sol.t)):
                time = integrated_sol.t[t_index]

                # calculate dq_qt for RC and d2q_dt2 for LCR
                q = dict()

                i = 0
                for wire in wires:
                    if wire.device.inductance is not None:
                        q[wire] = (integrated_sol.y[i][t_index], integrated_sol.y[i + 1][t_index])
                        i += 2
                    else:
                        q[wire] = (integrated_sol.y[i][t_index],)
                        i += 1

                derivatives = self.__derivatives(time, q, first_law, second_law)

                # add all values to solution
                for wire in wires:
                    time_list, q_list = solution[wire]
                    time_list.append(time)
                    q_list.append((*q[wire], derivatives[wire]))

                # current and its derivative is 0 for void_wires, charge is init_charge or last known value
                for wire in self.__voidWires:
                    time_list, q_list = solution[wire]
                    time_list.append(time)
                    if len(q_list) == 0:
                        charge = wire.device.initCharge
                    else:
                        charge = q_list[-1][0]

                    if wire.device.inductance is not None:
                        q_list.append((charge, 0.0, 0.0))
                    else:
                        q_list.append((charge, 0.0))

            # call event functions
            for func in event[1]:
                func[0](*func[1], **func[2])

            # recalculate effective wires and junctions
            self.__simplifyCircuit()
            degree = len(self.__effectiveWires)
            wires = random.sample(list(self.__effectiveWires), degree)

            # set init_conditions to the last calculated values
            init_conditions = []
            for wire in wires:
                time_list, q_list = solution[wire]
                init_conditions.extend(q_list[-1][:-1])

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

    def __get_intervals(self, end: float) -> List[Tuple[Tuple[float, float], List[Tuple[Callable, Any, Any]]]]:
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
