"""
<<<THIS FILE SHOWS AN EXAMPLE OF HOW TO FIND EFFECTIVE RESISTANCE OF A CIRCUIT>>>
"""
from viggy_circuits import *


def main():
    circuit = Circuit()

    # Wires of a wheatstone bridge

    # first pair of arms
    wireAB = Wire(resistance=1)
    wireAC = Wire(resistance=4)

    # second pair of arms
    wireDB = Wire(resistance=3)
    wireDC = Wire(resistance=12)

    # the galvanometer wire
    wireBC = Wire(resistance=4)

    junctionA = Junction()
    junctionB = Junction()
    junctionC = Junction()
    junctionD = Junction()

    # make the connections
    circuit.connect(junctionA, junctionB, wireAB)
    circuit.connect(junctionA, junctionC, wireAC)
    circuit.connect(junctionD, junctionB, wireDB)
    circuit.connect(junctionD, junctionC, wireDC)
    circuit.connect(junctionB, junctionC, wireBC)

    # a battery with some internal resistance is connected between A and D
    battery = Wire(resistance=1, battery=1)

    circuit.connect(junctionA, junctionD, battery)

    # find the currents
    solution = circuit.initialCurrents()

    i = solution[battery]
    V = battery.potentialDrop(i=i)

    print(f"Resistance of wheatstone bridge is {abs(V / i)}")


if __name__ == '__main__':
    main()
