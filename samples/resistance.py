"""
<<<THIS FILE SHOWS AN EXAMPLE OF HOW TO FIND EFFECTIVE RESISTANCE OF A CIRCUIT>>>
"""
from viggy_circuits import *

circuit = Circuit()

# Wires of a wheatstone bridge

# first pair of arms
wireAB = Wire(Device(resistance=1))
wireAC = Wire(Device(resistance=4))

# second pair of arms
wireDB = Wire(Device(resistance=3))
wireDC = Wire(Device(resistance=12))

# the galvanometer wire
wireBC = Wire(Device(resistance=4))

# the battery (with voltage 1)
V = 1
battery = Wire(Device(resistance=0.1,
                      battery=V))

junctionA = Junction()
junctionB = Junction()
junctionC = Junction()
junctionD = Junction()

# make the connections
circuit.connect(junctionA, battery, wireAB, wireAC)
circuit.connect(junctionD, battery, wireDB, wireDC)
circuit.connect(junctionB, wireAB, wireDB, wireBC)
circuit.connect(junctionC, wireAC, wireDC, wireBC)

solution = circuit.solve(end=1, dt=1)
