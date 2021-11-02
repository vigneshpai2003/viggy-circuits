from circuits import *
from matplotlib import pyplot as plt

# create the global circuit object
circuit = Circuit()

# create a list of wires
"""
each wire must be initialized using a Device

"""
wires = [Wire(Device(resistance=lambda t: 1.5 + sin(t),
                     battery=lambda t: 2.0,
                     capacitance=lambda t: 1.0,
                     inductance=lambda t: 0.2)),
         Wire(Device(resistance=lambda t: 0.1,
                     battery=AC(1.0, 5.0, 0.5),
                     capacitance=lambda t: 1.0)),
         Wire(Device(resistance=lambda t: 1.1,
                     capacitance=lambda t: 1.5 + sin(t),
                     init_charge=-5.0,
                     inductance=lambda t: 1.0)),
         Wire(Device(resistance=lambda t: 1.0,
                     capacitance=lambda t: 1.0,
                     inductance=lambda t: 1.5 + sin(t)))]
