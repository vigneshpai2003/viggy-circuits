# Circuits

This module simulates arbitrarily large electric circuits and provides the 
solution (i.e. current) in each wire of the circuit.

### Dependencies

``numpy``, ``scipy.linalg``, ``scipy.integrate``

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
`init_charge` and `init_current` respectively may also be given,
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
circuit.connect(junction, {wire1, wire2, wire3})
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

The solution is a dictionary with the keys being the wire,
the value is a tuple of two lists.

The first list is a list of time values at which the solution was evaluated.
The second list is a list of <img src="https://render.githubusercontent.com/render/math?math=(q, i, \frac{di}{dt})"> tuples,
note that in case inductance was not given for wire, it would be a list of <img src="https://render.githubusercontent.com/render/math?math=(q, i)"> tuples.

We can extract the values into individual values using the `zip` function:

```
time = solution[wire][0]
q, i, di_dt = list(zip(*solution[wire][1]))
```

### Visualization of Solution
```
from matplotlib import pyplot as plt

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
```
![Matplotlib plot](https://github.com/vigneshpai2003/Circuits/blob/master/plots/rlc1.png)
