from viggy_circuits import *
from fractions import Fraction


def R(n):
    circuit = Circuit()
    polyA = [Junction() for _ in range(n)]
    polyB = [Junction() for _ in range(n)]

    polyAWires = [Wire(Device(resistance=1)) for _ in range(n)]
    polyBWires = [Wire(Device(resistance=1)) for _ in range(n)]
    ABWires = [Wire(Device(resistance=1)) for _ in range(n)]

    # make the connections of the prism
    for j in range(n):
        circuit.connect(polyA[j], polyAWires[j], polyAWires[j-1], ABWires[j])
        circuit.connect(polyB[j], polyBWires[j], polyBWires[j-1], ABWires[j])

    battery = Wire(Device(resistance=1, battery=1))

    # connect battery to first two vertices of A
    circuit.connect(polyA[0], battery)
    circuit.connect(polyA[1], battery)

    i = circuit.initialCurrents()[battery]
    V = battery.device.potentialDrop(t=0, i=i)

    return abs(V / i)


if __name__ == '__main__':
    for i in range(3, 20):
        print(f"R({i}) = {Fraction(R(i)).limit_denominator()}")
