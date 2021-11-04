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
Wire objects may be created using Device which takes several arguments

+ `resistance`
+ `battery`
+ `capacitance`
+ `inductance`

These arguments may be floats or functions of time that return floats.
Resistance is mandatory while the rest are optional,
there must be a battery in at least one wire in the circuit.

In case capacitance or inductance is given,
`initCharge` and `initCurrent` respectively may also be given,
by default they are `0.0`.

Example:
```
wire = Wire(Device(resistance = 1.0,
                   capacitance = lambda t: 3.0 + math.sin(t)))
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
+ Each junction is connected to at least two wires
+ In case either of the above rules are violated,
  the simulation will ignore the violators and will not
  consider them to be part of the circuit.

### Solving and Analysis of Solution

The solution may be calculated till the end time with given time step:

```
solution = circuit.solve(end=20, dt=0.03)
```

The solution of a specific wire can be accessed as if solution is a dictionary,
the solution for a wire can be unpacked to give arrays of some quantities.

```
wireSolution = solution[wire]

t, q, i, di_dt = wireSolution
```

Note: if `wire` does not have inductor, there would not be `di_dt` in solution

### Visualization of Solution
```
from matplotlib import pyplot as plt

plt.figure(num="viggy-circuits")
plt.title("LCR Circuit")
plt.xlabel("time (s)")
plt.xticks(range(end + 1))

plt.plot(t, q)
plt.plot(t, i)
plt.plot(t, di_dt)

plt.legend(['q', 'i', 'di/dt'])
plt.grid(True)
plt.show()
```
![Matplotlib plot](https://github.com/vigneshpai2003/viggy-circuits/blob/master/plots/plot1.png?raw=True)
