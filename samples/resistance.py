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

junctionA = Junction()
junctionB = Junction()
junctionC = Junction()
junctionD = Junction()

# make the connections
circuit.connect(junctionA, wireAB, wireAC)
circuit.connect(junctionD, wireDB, wireDC)
circuit.connect(junctionB, wireAB, wireDB, wireBC)
circuit.connect(junctionC, wireAC, wireDC, wireBC)

# a battery with some internal resistance is connected between A and D
battery = Wire(Device(resistance=1, battery=1))

circuit.connect(junctionA, battery)
circuit.connect(junctionD, battery)

# find the currents
solution = circuit.initialCurrents()

i = solution[battery]
V = battery.device.potentialDrop(t=0, i=i)

print(f"Resistance of wheatstone bridge is {abs(V / i)}")
