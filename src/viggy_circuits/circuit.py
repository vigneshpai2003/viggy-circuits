from __future__ import annotations

import random

import numpy as np
from scipy import linalg, integrate
import networkx as nx
import matplotlib.pyplot as plt

from .wireCollection import WireCollection
from .solution import Solution

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Set, List, Dict, Tuple, Callable, Any

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

    def __firstLawEquations(self) -> Set[WireCollection]:
        """
        if there are n Junctions, the sum of the equations formed by them is 0 (linearly dependent)
        this is because each wire shows up twice but with different sign each time
        It can also be proved that any n-1 Junction equations are linearly independent
        :return: the wires to be used in implementation of Kirchhoff's First Law
        """
        equations: Set[WireCollection] = set()

        junctionIter = iter(self.__effectiveJunctions)
        next(junctionIter)  # we dont need one junction

        while True:
            try:
                junction = next(junctionIter)
                junctionWires = junction.wires - self.__voidWires

                """
                in case all are LCR wires then all will be used in writing first law equation with di_dt
                in case only some are LCR, we already know their i values and hence they are not variables
                """
                allLCR = all([wire.isLCR for wire in junctionWires])
                equations.add(WireCollection(*[(wire, wire.sign(junction))
                                               for wire in junctionWires if not wire.isLCR or allLCR]))

            except StopIteration:
                return equations

    def __secondLawEquations(self) -> Set[WireCollection]:
        """
        :return: the wires to be used in implementation of Kirchhoff's Second Law;
        """
        # create a multi graph of effective circuit
        graph = nx.MultiGraph()

        for wire in self.__effectiveWires:
            graph.add_edge(*wire.junctions, wire=wire)

        """
        cycle_basis algorithm only works on Graph and not MultiGraph
        for the junctions in between which more that one wires are present, we take some one, say X
        for the rest of the parallel wires, we can apply KVL between the parallel wire and X
        """
        pseudoGraph = nx.Graph(graph)

        equations = set()
        parsedParallelWires = set()

        for cycle in nx.algorithms.cycle_basis(nx.Graph(graph)):
            equation = WireCollection()  # equation corresponding to cycle
            for i in range(-1, len(cycle) - 1):
                junction1, junction2 = cycle[i], cycle[i + 1]
                wire = pseudoGraph.edges[junction1, junction2]["wire"]
                equation.append(wire, wire.sign(junction1))

                # get all the parallel wires that havent already been parsed
                parallelWires = set()

                for d in graph.get_edge_data(junction1, junction2).values():
                    parallelWire = d['wire']
                    if parallelWire is not wire and parallelWire not in parsedParallelWires:
                        parallelWires.add(parallelWire)

                # add the KVL equations for the parallelWires
                for parallelWire in parallelWires:
                    equations.add(WireCollection((wire, wire.sign(junction1)),
                                                 (parallelWire, parallelWire.sign(junction2))))

                # mark these parallel wires as parsed
                parsedParallelWires.update(parallelWires)

            equations.add(equation)

        return equations

    @staticmethod
    def __derivatives(t: float, x: np.ndarray, wires: List[Wire], wireToIndex: Dict[Wire, int],
                      firstLawEquations: Set[WireCollection],
                      secondLawEquations: Set[WireCollection]) -> np.ndarray:
        """
        :param t: time
        :param x: a list of charge and current values
        :param wires: the list of wires corresponding to each value in x
        :param wireToIndex: the index of the charge / current value in x for each wire

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
                    V[i] -= sign * x[wireToIndex[wire] + 1]
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
                    V[i] -= sign * wire.device.resistance(t) * x[wireToIndex[wire] + 1]
                else:
                    R[i, wires.index(wire)] = sign * wire.device.resistance(t)

                # potential drop due to battery
                if wire.device.battery is not None:
                    V[i] += sign * wire.device.battery(t)

                # potential drop due to capacitor
                if wire.device.capacitance is not None:
                    V[i] -= sign * x[wireToIndex[wire]] / wire.device.capacitance(t)

            i += 1

        try:
            derivatives = linalg.solve(R, V)
        except linalg.LinAlgError:
            print(firstLawEquations, secondLawEquations)

        dx_dt = np.zeros_like(x)

        for wire in wires:
            if wire.isLCR:
                dx_dt[wireToIndex[wire]] = x[wireToIndex[wire] + 1]  # i is already in input
                dx_dt[wireToIndex[wire] + 1] = derivatives[wires.index(wire)]  # di_dt
            else:
                dx_dt[wireToIndex[wire]] = derivatives[wires.index(wire)]  # i

        return dx_dt

    def solve(self, end: float, dt: float) -> Solution:
        solution = Solution(self.__wires)

        for (start, stop), events in self.__getIntervals(end):
            # calculate effective wires and junctions
            self.__simplifyCircuit()
            degree = len(self.__effectiveWires)
            wires = random.sample(list(self.__effectiveWires), degree)

            # set initial conditions of simulation
            initConditions = solution.lastKnownValues(wires)
            indexOf = solution.mapWireToIndex(wires)

            # calculate firstLaw and secondLaw
            firstLaw = self.__firstLawEquations()
            secondLaw = self.__secondLawEquations()

            # method: Nonstiff = ('RK45', 'RK23', 'DOP853'); Stiff = ('Radau', 'BDF'); Universal = 'LSODA'
            integral = integrate.solve_ivp(fun=Circuit.__derivatives, y0=initConditions,
                                           t_span=(start, stop), t_eval=np.arange(start, stop, dt),
                                           args=(wires, indexOf, firstLaw, secondLaw))

            # update solution
            for i in range(len(integral.t)):
                t = integral.t[i]
                x = integral.y[:, i]

                derivatives = self.__derivatives(t, x, wires, indexOf, firstLaw, secondLaw)

                # add all values to solution
                for wire in wires:
                    solution[wire].update(t, q=x[indexOf[wire]],
                                          i=derivatives[indexOf[wire]],
                                          di_dt=derivatives[indexOf[wire] + 1] if wire.isLCR else None)

                for wire in self.__voidWires:
                    solution[wire].updateVoid(t)

            # call event functions
            for func, args, kwargs in events:
                func(*args, **kwargs)

        solution.freeze()

        return solution

    def addEvent(self, time: float, event: Callable, args=None, kwargs=None):
        """
        when solving, the event is called using args and kwargs at specified time;
        multiple events are allowed at same time, but add_event must be called multiple times;
        """
        if time < 0:
            raise ValueError(f"time must be non-negative, cannot be {time}")

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
        splits time from 0 to end into smaller intervals separated by event execution
        :return: List[( (start, stop), List[(func, args, kwargs)] )]
        """
        times = sorted([i for i in list(self.__events.keys()) if i <= end])

        if len(times) == 0:
            return [((0, end), [])]

        intervals = []
        old = 0.0
        for i in range(len(times)):
            intervals.append(((old, times[i]), self.__events[times[i]]))
            old = times[i]

        intervals.append(((times[-1], end), []))

        return intervals

    def initialCurrents(self) -> Dict[Wire, float]:
        """
        :return: a dictionary mapping each wire to current flowing in it
        """
        # calculate effective wires and junctions
        self.__simplifyCircuit()
        degree = len(self.__effectiveWires)
        wires = random.sample(list(self.__effectiveWires), degree)

        x = Solution.initialValues(wires)
        indexOf = Solution.mapWireToIndex(wires)

        # calculate firstLaw and secondLaw
        firstLaw = self.__firstLawEquations()
        secondLaw = self.__secondLawEquations()

        dx_dt = self.__derivatives(0, x, wires, indexOf, firstLaw, secondLaw)

        return dict([(wire, dx_dt[indexOf[wire]]) for wire in wires])

    def showGraph(self):
        self.__simplifyCircuit()
        graph = nx.MultiGraph()

        for wire in self.__effectiveWires:
            graph.add_edge(*wire.junctions)

        nx.draw(nx.convert_node_labels_to_integers(graph), with_labels=True, font_weight='bold')
        plt.show()
