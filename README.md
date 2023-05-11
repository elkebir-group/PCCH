# Parsimonious Consistent Comigration History

Parsimonious Consistent Comigration History (PCCH) is a computational framework for inferring migration history from a given phylogeny with leaves labeled with locations.

# Table of Contents
1. [Dependencies](##dependencies)
2. [Usage instruction](##usage-instruction)

## Dependencies

### Gurobi

`PCCH` requires a valid Gurobi installation and license key. The location of Gurobi should be present in `LD_LIBRARY_PATH` (linux) and `DYLD_LIBRARY_PATH` (macOS) the license key should be saved in the environment variable `GRB_LICENSE_KEY`.

### Python

`PCCH` requires Python version 3.6 or newer.

### `networkx`

`PCCH` requires python package `networkx` to run. It is available with all python package managers.

        $ # install with pip
        $ pip install networkx[default]
        $
        $ # install with conda
        $ conda install -c anaconda networkx




## Usage Instruction

### I/O formats

`PCCH` is a program that takes as input two files: a tree file and a leaf labeling file. The tree file contains a list of edges, defining the structure of a tree, and the leaf labeling file contains the labels assigned to each edge. 

1. **Tree file** : The tree file contains a list of edges that define the structure of a tree. Each line in the file represents an edge, and the edges should be in the format: `vertex1 vertex2`. For example:

        1   2
        2   3
        2   4 
        1   3

2. **Leaf labeling file** : The leaf labeling file contains the labels assigned to each leaf. Each line in the file corresponds to a leaf of the tree represented in the tree file, and the labels should be in the format: `leaf label`. For example:

        1   A
        2   B
        3   C

`PCCH` outputs a location labeling file, which follows the same format as the leaf labeling, except it contains the labels of all the vertices of the tree represented in the tree file. It also returns a DOT file containing the tree along with location labeling. `PCCH` also prints `<primary location> <number of migrations> <number of comigrations> Optimal <running time (in seconds)>` on console.

### Usage

`PCCH` can be run using python.

        usage: PCCH.py [-h] [-p PRIMARY] [-c COLORMAP] [--log] [-o OUTPUT] [-t THREADS]
                    phylogeny leaf_labeling

        PCCH

        positional arguments:
        phylogeny             Input phylogeny
        leaf_labeling         Input leaf labeling

        optional arguments:
        -h, --help            show this help message and exit
        -p PRIMARY, --primary PRIMARY
                                Primary location
        -c COLORMAP, --colormap COLORMAP
                                Color map file
        --log                 Outputs Gurobi logging
        -o OUTPUT, --output OUTPUT
                                Output folder
        -t THREADS, --threads THREADS
                            Number of threads

An example execution

        $ python PCCH.py t.tree l.labeling
        a-      7       6       Optimal         1.366915225982666
        