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
Then, instantiate with the desired material, crystalline system, and orientation.
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

Example (found in examples.py). Note that all orientations are given as Python lists, while materials and crystal systems are given as Python strings:
```
from crystal import Crystal

#We first define a crystal named c1, composed of GaN in the cubic system (zincblende), of the class Crystal. 
#with Miller indices (001)[110]
c1 = Crystal(mat='GaN', system = 'cub', z=[0,0,1], x = [1,1,0])

#We add the 113 pole
c1.add_pole([1,1,3], None, None)
#and save it under the name "113 GaN"
c1.save_fig('113 GaN')

#We define another crystal named c2, composed of GaN in the cubic system (zincblende), of the class Crystal. 
#with Miller indices (001)[110]
c2 = Crystal(mat='GaN', system = 'cub', z=[0,0,1], x = [1,1,0])

#We add the four twins about the <110> axes by +- 109.5 degrees.
c2.add_pole(ref=[1,1,3], twin_axis=[1,-1,0], twin_angle=109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,-1,0], twin_angle=-109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,1,0], twin_angle=109.5)
c2.add_pole(ref=[1,1,3], twin_axis=[1,1,0], twin_angle=-109.5)

#and save the file under the name "Twinned GaN".
c2.save_fig('Twinned GaN')

```


