from __future__ import annotations

import copy
from random import sample

from numpy import arange
from scipy import linalg, integrate

from .error import CircuitError
from .wireCollection import WireCollection, Direction

from typing import TYPE_CHECKING

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
        sorts the wires into effective wires and null wires (which do not contribute to calculation);
        similarly for junctions
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

    def __get_loops(self, base_wire: Wire) -> List[WireCollection]:
        """
        :param base_wire: should be in self.__effectiveWires
        :return: all the closed loops with the base_wire;
        """
        loops = []
        endpoint = base_wire.junction1

        def __loops(current_junction: Junction, current_wire: Wire,
                    loop_path: WireCollection, parsed_junctions: list):
            """
            recursive function that calculates closed loops
            """
            for new_wire in current_junction.otherWires(current_wire) - self.__voidWires:
                """
                for every wire in the junction,
                an attempt is made to add it to loop_path
                if failed, the entire loop path corresponding up to that point is discarded
                """
                if new_wire.isNull:
                    continue

                # direction is assumed to be from junction1 to junction2 always
                if current_junction is new_wire.junction1:
                    sign = Direction.FORWARD
                elif current_junction is new_wire.junction2:
                    sign = Direction.BACKWARD
                else:
                    raise CircuitError("sign could not be deduced")

                # new_loop_path must be copied and not referenced
                new_loop_path = copy.copy(loop_path)
                new_loop_path.append(new_wire, sign)

                next_junction = new_wire.otherJunction(current_junction)

                if next_junction is endpoint:
                    # if closed loop is found
                    loops.append(new_loop_path)
                    continue
                elif next_junction in parsed_junctions:
                    # if junction was already encountered before
                    continue

                __loops(current_junction=next_junction, current_wire=new_wire,
                        loop_path=new_loop_path, parsed_junctions=[*parsed_junctions, next_junction])

        __loops(current_junction=base_wire.junction2, current_wire=base_wire,
                loop_path=WireCollection((base_wire, Direction.FORWARD)), parsed_junctions=[base_wire.junction2])

        return loops

    def __first_law_wires(self) -> Set[WireCollection]:
        """
        return the wires to be used in implementation of Kirchhoff's First Law;
        """
        equations = set()
        # possible_equations is updated to include all possible combinations of already formed equations
        possible_equations = set()

        for junction in self.__effectiveJunctions:
            # remove wires that are void
            junction_wires = junction.wires - self.__voidWires

            if all([wire.device.inductance is not None for wire in junction_wires]):
                # if all wires have inductor
                effective_junction_wires = junction_wires.copy()
            else:
                # if some or no wires have inductor
                # current in all wires with inductors are known in __derivatives
                # only those without inductor should be included
                effective_junction_wires = set()
                for wire in junction_wires:
                    if wire.device.inductance is None:
                        effective_junction_wires.add(wire)

            # check if these wire can be derived from existing wire equations
            if effective_junction_wires in possible_equations:
                continue

            # update possible combos
            possible_equations.add(frozenset(effective_junction_wires))
            possible_equations.update(set(wire_combo.symmetric_difference(effective_junction_wires)
                                          for wire_combo in possible_equations))
            """
            if 2 equation in the currents of a certain set of wire exists,
            another equation can be created by subtracting them,
            this new equation will not contain variables common to both equations,
            it will contain the symmetric difference of the variables
            Note: the orientation of wire is not being considered, should it be?
            """

            loop = WireCollection()
            for wire in junction_wires:
                if junction is wire.junction1:
                    loop.append(wire, Direction.FORWARD)
                else:
                    loop.append(wire, Direction.BACKWARD)
            equations.add(loop)

        return equations

    def __second_law_wires(self, limit: int) -> Set[WireCollection]:
        """
        return the wires to be used in implementation of Kirchhoff's Second Law;
        """
        equations = set()
        parsed_wires = set()

        for wire in self.__effectiveWires:
            for loop in self.__get_loops(base_wire=wire):
                if set(loop.wires).issubset(parsed_wires):
                    # if equation in these variablles is already formed
                    continue
                parsed_wires.update(set(loop.wires))
                equations.add(loop)
                if len(equations) == limit:
                    return equations

        return equations

    def __derivatives(self, t: float, q: Dict[Wire, Union[Tuple[float], Tuple[float, float]]],
                      first_law_equations: Set[WireCollection],
                      second_law_equations: Set[WireCollection]) -> Dict[Wire, float]:
        """
        t is time and x is dictionary corresponding each wire to their respective charges

        q is of the form [q,] for wires without inductor
        q is of the form [q, dq_qt] for wires with inductor

        returns dq_dt for wires without inductor and d2q_dt2 for wires with inductor in the form of a dictionary
        """
        degree = len(self.__effectiveWires)  # number of variables to solve for
        wires = sample(tuple(self.__effectiveWires), degree)  # this will indicate order of wires in matrix

        # R i = V
        R = []
        V = []

        # fill R, V using first law
        for equation in first_law_equations:
            R_eq = [0.0] * degree
            V_eq = 0.0

            if all([wire.device.inductance is not None for wire in equation.wires]):
                # all are LCR
                for wire in equation.wires:
                    R_eq[wires.index(wire)] = float(equation.getDirection(wire))
            else:
                # if atleast one is RC
                for wire in equation.wires:
                    if wire.device.inductance is not None:
                        # if LCR then dq_dt is known
                        V_eq -= float(equation.getDirection(wire)) * q[wire][1]
                    else:
                        # if RC
                        R_eq[wires.index(wire)] = float(equation.getDirection(wire))

            R.append(R_eq)
            V.append(V_eq)

        # fill R, V using second law
        for equation in second_law_equations:
            R_eq = [0.0] * degree
            V_eq = 0.0

            for wire in equation.wires:
                sign = float(equation.getDirection(wire))

                if wire.device.inductance is not None:
                    # if LCR
                    R_eq[wires.index(wire)] = sign * wire.device.inductance(t)

                    # calculate potential drop due to resistor
                    V_eq -= sign * wire.device.resistance(t) * q[wire][1]
                else:
                    # if RC
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
        wires = sample(list(self.__effectiveWires), degree)

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
            first_law = self.__first_law_wires()
            second_law = self.__second_law_wires(limit=degree - len(first_law))

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
            wires = sample(list(self.__effectiveWires), degree)

            # set init_conditions to the last calculated values
            init_conditions = []
            for wire in wires:
                time_list, q_list = solution[wire]
                init_conditions.extend(q_list[-1][:-1])

        return solution

    def add_event(self, time: float, event: Callable, args=None, kwargs=None):
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
