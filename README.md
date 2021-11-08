# viggy-circuits

This module simulates arbitrarily large electric circuits and provides the 
solution (i.e. current) in each wire of the circuit.

### Installation

`> git clone https://github.com/vigneshpai2003/viggy-circuits.git`

`> cd viggy-circuits`

`> pip install -r requirements.txt`

`> pip install .`

Try running any of the example in samples folder

### Basic Tutorial

The Circuit object will be responsible for all calculations

```
circuit = Circuit()
```

#### Wire Object
Wire objects take several arguments

+ `resistance`
+ `battery`
+ `capacitance`
+ `inductance`

These arguments may be floats or functions of time that return floats.
Resistance is mandatory for ResistorWire.
Inductance is mandatory for InductorWire.

In case capacitance or inductance is given,
`initCharge` and `initCurrent` respectively may also be given,
by default they are `0.0`.

Example:
```
wire = ResistorWire(resistance = 1.0,
                    capacitance = lambda t: 3.0 + math.sin(t))
```

#### Connecting Wires
Wires are connected via Junction objects which may be initialized as:

```
junction = Junction()
```

The connection is made using the circuit object

```
circuit.connect(junction, wire1, wire2, wire3)
```

Note:

+ Each wire is connected to exactly two junctions
+ Each junction is connected to several wires
+ Note that if a junction has 0 or 1 wires, no current will flow through it

### Solving and Analysis of Solution

The solution may be calculated till the end time with given time step,
the solution of a specific wire can be accessed as if solution is a dictionary

```
solution = circuit.solve(end=20, dt=0.01)
wireSolution = solution[wire]
```

### Visualization of Solution
```
from matplotlib import pyplot as plt

plt.figure(num="viggy-circuits")
plt.title("LCR Circuit")
plt.xlabel("time (s)")
plt.xticks(range(end + 1))

plt.plot(wireSolution.t, wireSolution.q)  # charge
plt.plot(wireSolution.t, wireSolution.i)  # current
plt.plot(wireSolution.t, wireSolution.di_dt)  # current gradient
plt.plot(wireSolution.t, wireSolution.v)  # potential drop
plt.plot(wireSolution.t, wireSolution.powerR)  # power dissipated across Resistor
plt.plot(wireSolution.t, wireSolution.energyC)  # energy in capacitor
plt.plot(wireSolution.t, wireSolution.energyL)  # energy in inductor

plt.legend(['$q$', '$i$', r'$\frac{di}{dt}$', '$V$', '$P_R$', '$U_C$', '$U_L$'])
plt.grid(True)
plt.show()

```
![plot](https://github.com/vigneshpai2003/viggy-circuits/blob/master/plots/plot2.png?raw=True)
