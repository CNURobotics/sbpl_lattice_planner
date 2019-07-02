
Python-based Motion Primitive Generators
======

This folder contains a collection of simple Python 2.7-based scripts for generating and visualizing motion primitives.  These Python scripts are based on Matlab scripts released by Maxim Likhachev as part of the SBPL library (Copyright (c) 2008).

### Main scripts

#### `mprim_generator.py`

This script generates the motion primitive file from a simplified input file.
See `mprim_generator_050.csv` and `mprim_generator_100.csv` for example input files.

* Sample Usage: `python2 mprim_generator.py -i mprim_generator_050.csv -o mprim_050.prim -r 0.05 -a`

This generates the motion primitive file including some Arc-Line-Arc swerve motions.

Use the `--help` option to see all arguments.

The primitive definitions work with the SBPL library that allows for non-uniform angles that align with the grid cells, and a custom  branch with a modification that specifies the action cost as part of the primitive.


#### `mprim_grid_gen.py`

This script generates a list of potential primitive targets by evaluating a grid of data to see which ones allow for valid primitives.

Typical usage is to pipe the output to a file and use to manually select nodes to define the generator files (e.g. `mprim_generator_050.csv`)

This script is less developed and has more hard coded values, but is useful.

### `mprim_visualize.py`

This file loads a primitive definition file and visualizes the motion primitives using MatPlotLib.

 Sample Usage: `python2 mprim_visualize.py -i mprim_050.prim`


### Utility files

The following files contain functions used by the other scripts:
* `mprim_generator_generators.py`
* `mprim_generator_io.py`
* `mprim_generator_utils.py`
