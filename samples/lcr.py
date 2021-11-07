"""
<<<THIS FILE SHOWS AN EXAMPLE OF HOW TO SIMULATE AN LCR CIRCUIT>>>
"""
from viggy_circuits import *

from math import sin
import matplotlib.pyplot as plt


def main():
    circuit = Circuit()

    wires = [Wire(Device(resistance=lambda t: 1.5 + sin(t),
                         battery=2.0,
                         capacitance=1.0,
                         inductance=0.2)),
             Wire(Device(resistance=0.1,
                         battery=generatorAC(emf=1.0, f=1.0, phi=0.5),  # AC generates an alternating emf function
                         capacitance=1.0)),
             Wire(Device(resistance=1.1,
                         capacitance=lambda t: 1.5 + sin(t),
                         initCharge=-5.0,
                         inductance=1.0)),
             Wire(Device(resistance=1.0,
                         capacitance=1.0,
                         inductance=lambda t: 1.5 + sin(t)))]

    junctions = [Junction(), Junction(), Junction()]

    circuit.connect(junctions[0], wires[0], wires[2], wires[3])
    circuit.connect(junctions[1], wires[0], wires[1])
    circuit.connect(junctions[2], wires[1], wires[2], wires[3])

    def toggle(_index):
        # function that toggles the switch of the _index wire
        wires[_index].switch.toggle()

    # events are functions that are executed with specified arguments at specified time in the simulation
    circuit.addEvent(time=3.0, event=toggle, args=(3,))
    circuit.addEvent(time=7.0, event=toggle, args=(2,))
    circuit.addEvent(time=7.0, event=toggle, args=(3,))
    circuit.addEvent(time=8.0, event=toggle, args=(2,))

    # circuit.showGraph()

    end = 20
    solution = circuit.solve(end=end, dt=0.01)
    wireSolution = solution[wires[2]]

    # plot solution using matplotlib
    # plt.figure(num="viggy-circuits")
    # plt.title("LCR Circuit")
    # plt.xlabel("time (s)")
    # plt.xticks(range(end + 1))
    #
    # plt.plot(wireSolution.t, wireSolution.q)  # charge
    # plt.plot(wireSolution.t, wireSolution.i)  # current
    # plt.plot(wireSolution.t, wireSolution.di_dt)  # current gradient
    # plt.plot(wireSolution.t, wireSolution.v)  # potential drop
    # plt.plot(wireSolution.t, wireSolution.powerR)  # power dissipated across Resistor
    # plt.plot(wireSolution.t, wireSolution.energyC)  # energy in capacitor
    # plt.plot(wireSolution.t, wireSolution.energyL)  # energy in inductor
    #
    # plt.legend(['$q$', '$i$', r'$\frac{di}{dt}$', '$V$', '$P_R$', '$U_C$', '$U_L$'])
    # plt.grid(True)
    # plt.show()


if __name__ == '__main__':
    while True:
        main()
        print("done")
