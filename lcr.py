"""
<<<THIS FILE SHOWS AN EXAMPLE OF HOW TO SIMULATE AN LCR CIRCUIT>>>
"""
from circuits import *
from matplotlib import pyplot as plt

circuit = Circuit()

wires = [Wire(Device(resistance=lambda t: 1.5 + sin(t),
                     battery=2.0,
                     capacitance=1.0,
                     inductance=0.2)),
         Wire(Device(resistance=0.1,
                     battery=AC(emf_rms=1.0, omega=5.0, phi=0.5),  # AC generates an alternating emf function
                     capacitance=1.0)),
         Wire(Device(resistance=1.1,
                     capacitance=lambda t: 1.5 + sin(t),
                     init_charge=-5.0,
                     inductance=1.0)),
         Wire(Device(resistance=1.0,
                     capacitance=1.0,
                     inductance=lambda t: 1.5 + sin(t)))]

junctions = [Junction(), Junction(), Junction()]

circuit.connect(junctions[0], {wires[0], wires[2], wires[3]})
circuit.connect(junctions[1], {wires[0], wires[1]})
circuit.connect(junctions[2], {wires[1], wires[2], wires[3]})


def toggle(_index):
    # function that toggles the switch of the _index wire
    wires[_index].switch.toggle()


# events are functions that are executed with specified arguments at specified time in the simulation
circuit.add_event(time=3.0, event=toggle, args=(3,))
circuit.add_event(time=7.0, event=toggle, args=(2,))
circuit.add_event(time=7.0, event=toggle, args=(3,))
circuit.add_event(time=8.0, event=toggle, args=(2,))

end = 20
solution = circuit.solve(end=end, dt=0.01)

time = solution[wires[2]][0]
q, i, di_dt = list(zip(*solution[wires[2]][1]))

# plot solution using matplotlib
plt.figure(num="Electric Circuits")
plt.title("LCR Circuit")
plt.xlabel("time (s)")
plt.xticks(range(end + 1))

plt.plot(time, q)
plt.plot(time, i)
plt.plot(time, di_dt)

plt.legend(['q', 'i', 'di/dt'])
plt.grid(True)
plt.show()
