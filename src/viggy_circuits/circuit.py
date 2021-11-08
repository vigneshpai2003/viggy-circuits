from __future__ import annotations

from collections import OrderedDict

import numpy as np
from scipy import linalg, integrate
import networkx as nx
import matplotlib.pyplot as plt

from .solution import Solution

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Set, List, Dict, Tuple, Callable, Any

    from .wire import Wire, Direction
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

    def connect(self, junction1: Junction, junction2: Junction, wire: Wire):
        self.__wires.add(wire)
        self.__junctions.update({junction1, junction2})

        junction1.connect(wire)
        junction2.connect(wire)

        wire.connect(junction1, junction2)

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

    def __firstLawEquations(self) -> List[List[Tuple[Wire, Direction]]]:
        """
        if there are n Junctions, the sum of the equations formed by them is 0 (linearly dependent)
        this is because each wire shows up twice but with different sign each time
        It can also be proved that any n-1 Junction equations are linearly independent
        :return: the wires to be used in implementation of Kirchhoff's First Law
        """
        equations: List[List[Tuple[Wire, Direction]]] = []

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
                equations.append([(wire, wire.direction(junction)) for wire in junctionWires if not wire.isLCR or allLCR])

            except StopIteration:
                return equations

    def __secondLawEquations(self) -> List[List[Tuple[Wire, Direction]]]:
        """
        :return: the wires to be used in implementation of Kirchhoff's Second Law;
        """
        # create a new sub effectiveGraph
        effectiveGraph = nx.MultiGraph()

        for wire in self.__effectiveWires:
            effectiveGraph.add_edge(*wire.junctions, wire=wire)

        """
        cycle_basis algorithm only works on Graph and not MultiGraph
        for the junctions in between which more that one wires are present, we take some one, say X
        for the rest of the parallel wires, we can apply KVL between the parallel wire and X
        """
        # converting multigraph to graph
        pseudoGraph = nx.Graph(effectiveGraph)

        equations: List[List[Tuple[Wire, Direction]]] = []
        parsedParallelWires: Set[Wire] = set()

        for cycle in nx.algorithms.cycle_basis(pseudoGraph):
            equation = []  # equation corresponding to cycle
            for i in range(len(cycle)):
                junction1, junction2 = cycle[i - 1], cycle[i]
                wire = pseudoGraph.edges[junction1, junction2]["wire"]
                equation.append((wire, wire.direction(junction1)))

                # get all the parallel wires that havent already been parsed
                parallelWires = set()

                for d in effectiveGraph.get_edge_data(junction1, junction2).values():
                    parallelWire = d['wire']
                    if parallelWire is not wire and parallelWire not in parsedParallelWires:
                        parallelWires.add(parallelWire)

                # add the KVL equations for the parallelWires
                for parallelWire in parallelWires:
                    equations.append([(wire, wire.direction(junction1)),
                                      (parallelWire, parallelWire.direction(junction2))])

                # mark these parallel wires as parsed
                parsedParallelWires.update(parallelWires)

            equations.append(equation)

        return equations

    @staticmethod
    def __derivatives(t: float, x: np.ndarray, wires: OrderedDict, wireToIndex: Dict[Wire, int],
                      firstLawEquations: List[List[Tuple[Wire, Direction]]],
                      secondLawEquations: List[List[Tuple[Wire, Direction]]]) -> np.ndarray:
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
            allLCR = all([wire.isLCR for wire, sign in equation])
            for wire, sign in equation:
                if wire.isLCR and not allLCR:
                    V[i] -= sign * x[wireToIndex[wire] + 1]
                else:
                    R[i, wires[wire]] = sign

            i += 1

        # fill R, V using second law
        for equation in secondLawEquations:
            for wire, sign in equation:
                if wire.isLCR:
                    R[i, wires[wire]] = sign * wire.inductance(t)

                    # calculate potential drop due to resistor
                    V[i] -= sign * wire.resistance(t) * x[wireToIndex[wire] + 1]
                else:
                    R[i, wires[wire]] = sign * wire.resistance(t)

                # potential drop due to battery
                if wire.battery is not None:
                    V[i] += sign * wire.battery(t)

                # potential drop due to capacitor
                if wire.capacitance is not None:
                    V[i] -= sign * x[wireToIndex[wire]] / wire.capacitance(t)

            i += 1

        derivatives = linalg.solve(R, V)

        dx_dt = np.zeros_like(x)

        for wire in wires:
            if wire.isLCR:
                dx_dt[wireToIndex[wire]] = x[wireToIndex[wire] + 1]  # i is already in input
                dx_dt[wireToIndex[wire] + 1] = derivatives[wires[wire]]  # di_dt
            else:
                dx_dt[wireToIndex[wire]] = derivatives[wires[wire]]  # i

        return dx_dt

    def solve(self, end: float, dt: float) -> Solution:
        solution = Solution(self.__wires)

        for (start, stop), events in self.__getIntervals(end):
            # calculate effective wires and junctions
            self.__simplifyCircuit()
            wires = OrderedDict([(wire, i) for i, wire in enumerate(self.__effectiveWires)])

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
        wires = OrderedDict([(wire, i) for i, wire in enumerate(self.__effectiveWires)])

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
