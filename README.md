# Pole_figure_simulation
Simulation of pole figure angles from a given material for any desired reflection. Rotations representing twins can be added.

Tested with:

 * Python 3.11.3 installed from https://www.python.org/ in Windows 11.
 * Python 3.8.10 included with the IDE Spyder installed from https://www.spyder-ide.org/ in Windows 11.

The following Python packages are needed: numpy, matplotlib, and pandas.
These packages can be installed through pip in a virtual environment (recommended) or global installation of Python.

crystal.py contains the main class powering the simulations. Two examples are given in examples.py for the calculation of 113 pole figure of zincblende GaN and
its twins around the <110> axes by +- 109.5 degress.

The simulations are stored in the folder SavedFigures. If it does not exist, the script will try to create it.

If more materials and crystal systems are desired, they can be added in database.py, which requires the six lattice parametters, a, b, c, alpha, beta, and gama.

The scripts are still under development.


# Usage:

In a new Python script, import the Crystal class from the crystal file:
```
from crystal import Crystal
```
Then, instantiate with the desired material, crystalline phase ('hex', 'cub', 'mon', <a>), and orientation.
```
c = Crystal(material, system, orientation_along_z, orientation_along_x)
```
Finally, add poles to the figure
```
c.add_pole(reflexion, twin_axis, twin_angle)
```
with twin_axis and twin_angle the information for twinning. Giving None to both will calculate the reflection for the orientation as is.

Finally, save the file:
```
c.save_fig(name)
```
wich will generate a png file in the folder SavedFigures.
