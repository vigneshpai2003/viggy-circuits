from __future__ import annotations

from typing import *
from enum import IntEnum
from math import sin
from random import sample

from scipy import linalg, integrate
from numpy import arange


class CircuitError(Exception):
    pass


class Switch:
    def __init__(self):
        self.__closed = True

    @property
    def is_closed(self) -> bool:
        return self.__closed

    @property
    def is_open(self) -> bool:
        return not self.is_closed

    def open(self):
        self.__closed = False

    def close(self):
        self.__closed = True

    def toggle(self):
        self.__closed = not self.__closed


class Device(object):
    def __init__(self, resistance: Union[Callable[[float], float], float],
                 battery: Union[Callable[[float], float], float] = None,
                 capacitance: Union[Callable[[float], float], float] = None,
                 inductance: Union[Callable[[float], float], float] = None,
                 init_charge=0.0, init_current=0.0):
        """
        :param resistance: a float or function of time that returns float
        :param battery: optional
        :param capacitance: optional
        :param inductance: optional
        :param init_charge: optional, only need to change if capacitance is given
        :param init_current: optional, only need to change if inductance is given
        """
        self.resistance = Device.process(resistance)
        self.battery = Device.process(battery)  # internal battery
        self.capacitance = Device.process(capacitance)
        self.inductance = Device.process(inductance)
        self.init_charge = init_charge
        self.init_current = init_current  # only significant if inductance is not None

    @classmethod
    def process(cls, x: Union[Callable[[float], float], float, None]) -> Union[Callable[[float], float], None]:
        if callable(x) or x is None:
            return x
        else:
            return lambda t: x

    def potential_drop(self, t: float, q: float = 0.0, dq_dt: float = 0.0, d2q_dt2: float = 0.0):

        v = 0.0

        if self.battery is not None:
            v += self.battery(t)

        if self.resistance is not None:
            v -= dq_dt * self.resistance(t)

        if self.capacitance is not None:
            v -= q / self.capacitance(t)

        if self.inductance is not None:
            v -= d2q_dt2 * self.inductance(t)

        return v


def AC(emf_rms: float, omega: float, phi: float) -> Callable[[float], float]:
    """
    :param emf_rms:
    :param omega:
    :param phi: phase constant
    :return: function of time that returns emf
    """

    def emf(t: float) -> float:
        return 2 ** 0.5 * emf_rms * sin(omega * t + phi)

    return emf


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

    def shallow_copy(self) -> WireCollection:
        new_obj = WireCollection()
        new_obj.__wires = [wire for wire in self.__wires]
        new_obj.__directions = [direction for direction in self.__directions]
        return new_obj

    def append(self, wire: Wire, direction: Direction):
        if wire in self.__wires:
            raise CircuitError(f"Wire: {wire} already exists in SignedWires: {self}")
        self.__wires.append(wire)
        self.__directions.append(direction)

    def get_direction(self, wire: Wire) -> Direction:
        return self.__directions[self.__wires.index(wire)]


class Circuit:
    def __init__(self):
        self.__wires: Set[Wire] = set()
        self.__effective_wires: Set[Wire] = set()
        self.__void_wires: Set[Wire] = set()
        self.__junctions: Set[Junction] = set()
        self.__effective_junctions: Set[Junction] = set()
        self.__void_junctions: Set[Junction] = set()
        self.__events = {}

    def connect(self, junction: Junction, wires: Set[Wire]):
        """
        function that connects junction to wire and each wire to junction
        """
        self.__junctions.add(junction)
        self.__wires.update(wires)

        junction.connect(*wires)
        for wire in wires:
            wire.connect(junction)

    def __simplify(self):
        """
        sorts the wires into effective wires and null wires (which do not contribute to calculation);
        similarly for junctions
        """
        self.__void_wires = set()
        self.__void_wires = set()

        for wire in self.__wires:
            if wire.is_null:
                self.__void_wires.add(wire)

        void_junction_count = -1
        while void_junction_count != 0:  # if no void junctions are found
            void_junction_count = 0
            # find void junctions
            for junction in self.__junctions - self.__void_junctions:
                if len(junction.wires - self.__void_wires) in (0, 1):
                    void_junction_count += 1
                    self.__void_junctions.add(junction)
                    self.__void_wires.update(junction.wires)

        self.__effective_wires = self.__wires - self.__void_wires
        self.__effective_junctions = self.__junctions - self.__void_junctions

    def __get_loops(self, base_wire: Wire) -> List[WireCollection]:
        """
        returns all the closed loops with the base_wire;
        base_wire should be an effective wire
        """
        loops = []
        endpoint = base_wire.junction1

        def __loops(current_junction: Junction, current_wire: Wire,
                    loop_path: WireCollection, parsed_junctions: list):
            """
            recursive function that calculates closed loops;
            """
            for new_wire in current_junction.other_wires(current_wire) - self.__void_wires:
                """
                for every wire in the junction,
                an attempt is made to add it to loop_path
                if failed, the entire loop path corresponding up to that point is discarded
                """
                if new_wire.is_null:
                    continue

                # direction is assumed to be from junction1 to junction2 always
                if current_junction is new_wire.junction1:
                    sign = Direction.FORWARD
                elif current_junction is new_wire.junction2:
                    sign = Direction.BACKWARD
                else:
                    raise CircuitError("sign could not be deduced")

                # new_loop_path must be copied and not referenced
                new_loop_path = loop_path.shallow_copy()
                new_loop_path.append(new_wire, sign)

                next_junction = new_wire.other_junction(current_junction)

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

        for junction in self.__effective_junctions:
            # remove wires that are void
            junction_wires = junction.wires - self.__void_wires

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

        for wire in self.__effective_wires:
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
        degree = len(self.__effective_wires)  # number of variables to solve for
        wires = sample(tuple(self.__effective_wires), degree)  # this will indicate order of wires in matrix

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
                    R_eq[wires.index(wire)] = float(equation.get_direction(wire))
            else:
                # if atleast one is RC
                for wire in equation.wires:
                    if wire.device.inductance is not None:
                        # if LCR then dq_dt is known
                        V_eq -= float(equation.get_direction(wire)) * q[wire][1]
                    else:
                        # if RC
                        R_eq[wires.index(wire)] = float(equation.get_direction(wire))

            R.append(R_eq)
            V.append(V_eq)

        # fill R, V using second law
        for equation in second_law_equations:
            R_eq = [0.0] * degree
            V_eq = 0.0

            for wire in equation.wires:
                sign = float(equation.get_direction(wire))

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
        self.__simplify()
        degree = len(self.__effective_wires)
        wires = sample(list(self.__effective_wires), degree)

        # set initial conditions
        init_conditions = []
        for wire in wires:
            if wire.device.inductance is not None:
                # if LCR
                init_conditions.extend([wire.device.init_charge, wire.device.init_current])
            else:
                # if RC
                init_conditions.append(wire.device.init_charge)

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
                for wire in self.__void_wires:
                    time_list, q_list = solution[wire]
                    time_list.append(time)
                    if len(q_list) == 0:
                        charge = wire.device.init_charge
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
            self.__simplify()
            degree = len(self.__effective_wires)
            wires = sample(list(self.__effective_wires), degree)

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
