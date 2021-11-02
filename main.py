from circuits import *
from matplotlib import pyplot


def main1():
    circuit = Circuit()
    wires = [Wire(device=Device(resistance=lambda t: 1.5 + sin(t),
                                battery=lambda t: 2.0,
                                capacitance=lambda t: 1.0,
                                inductance=lambda t: 0.2)),
             Wire(device=Device(resistance=lambda t: 0.1,
                                battery=AC(1.0, 5.0, 0.5),
                                capacitance=lambda t: 1.0)),
             Wire(device=Device(resistance=lambda t: 1.1,
                                capacitance=lambda t: 1.5 + sin(t),
                                init_charge=-5.0,
                                inductance=lambda t: 1.0)),
             Wire(device=Device(resistance=lambda t: 1.0,
                                capacitance=lambda t: 1.0,
                                inductance=lambda t: 1.5 + sin(t)))]

    junctions = [Junction(), Junction(), Junction()]

    circuit.connect(junctions[0], {wires[0], wires[2], wires[3]})
    circuit.connect(junctions[1], {wires[0], wires[1]})
    circuit.connect(junctions[2], {wires[1], wires[2], wires[3]})

    wires[0].swap_junctions()

    def toggle(_wire): wires[_wire].switch.toggle()

    circuit.add_event(3.0, event=toggle, args=(3,))
    circuit.add_event(7.0, event=toggle, args=(2,))
    circuit.add_event(7.0, event=toggle, args=(3,))
    circuit.add_event(8.0, event=toggle, args=(2,))

    end = 20
    solution = circuit.solve(end, 0.01)

    pyplot.figure(num="Electric Circuits")
    pyplot.title("LCR Circuit")
    pyplot.xlabel("time (s)")
    pyplot.xticks(range(end + 1))

    time = solution[wires[2]][0]
    q, dq_dt, d2q_dt2 = list(zip(*solution[wires[2]][1]))

    pyplot.plot(time, q)
    pyplot.plot(time, dq_dt)
    pyplot.plot(time, d2q_dt2)

    pyplot.legend(['Q', 'dQ/dt', 'd2Q/dt2'])
    pyplot.grid(True)
    pyplot.show()


def main2():
    pyplot.figure(num="Electric Circuits")
    pyplot.title("LCR Circuit")
    pyplot.xlabel("time (s)")
    pyplot.xticks(range(21))

    for i in range(5):

        circuit = Circuit()
        wires = [Wire(Device(resistance=lambda t: .2,
                             capacitance=lambda t: 1.0,
                             inductance=lambda t: 1.0,
                             init_charge=5 * i)),
                 Wire(Device(resistance=lambda t: 0.01,
                             battery=AC(0, 1, 0)))]

        junctions = [Junction(), Junction()]

        circuit.connect(junctions[0], set(wires))
        circuit.connect(junctions[1], set(wires))

        end = 20
        solution = circuit.solve(end, 0.01)

        time = solution[wires[0]][0]
        q, dq_dt, d2q_dt2 = list(zip(*solution[wires[0]][1]))

        pyplot.plot(time, dq_dt)

    pyplot.legend(['0', '1', '2', '3', '4'])
    pyplot.grid(True)
    pyplot.show()


if __name__ == '__main__':
    main2()
