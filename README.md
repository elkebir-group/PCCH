# TCC, PCC, and PCCH

The repository contains the tools, TCC, PCC, and PCCH, presented in the paper "Inferring Temporally Consistent Migration Histories" by  Mrinmoy Saha Roddur, Sagi Snir, and Mohammed El-Kebir. 

1. **Temporally Consistent Comigrations (TCC)**

TCC is a computational tool to check if a given location labeled phylogeny with the comigrations is temporally consistent or not. In our implementation, the tool `TCC.py` also incorporates `greedyComigrations` that infers comigrations from the given location labeled phylogeny. `TCC.py` takes as input a tree with the corresponding location labeling, runs `greedyComigrations` first to generate the comigrations, and then run TCC algorithm to check if the generated comigrations are temporally consistent or not. 

2. **Parsimonious Consistent Comigrations (PCC)**

Parsimonious Consistent Comigrations (PCC) is a computational tool for inferring temporally consistent comigrations from a given location labeled phylogeny.

3. **Parsimonious Consistent Comigration History**

Parsimonious Consistent Comigration History (PCCH) is a computational tool for inferring migration history from a given phylogeny with leaves labeled with locations.

# Table of Contents
1. [Installation](##installation)
2. [Usage instruction](##usage-instruction)

## Installation

### Python

`TCC`, `PCC`, and `PCCH` require Python version 3.6 or newer.

### `networkx`

`TCC`, `PCC`, and `PCCH` require python package `networkx` to run. It is available with all python package managers.

        $ # install with pip
        $ pip install networkx[default]
        $
        $ # install with conda
        $ conda install -c anaconda networkx

### Gurobi

`PCC` and `PCCH` require a valid Gurobi installation and license key. The location of Gurobi should be present in `LD_LIBRARY_PATH` (linux) and `DYLD_LIBRARY_PATH` (macOS) the license key should be saved in the environment variable `GRB_LICENSE_KEY`.

## Usage Instruction

### I/O formats

All of the tools take as input two files: a tree file and a location/leaf labeling file. The tree file contains a list of edges, defining the structure of a tree, and the location/leaf labeling file contains the labels assigned to each vertex/leaf. 

1. **Tree file** : The tree file contains a list of edges that define the structure of a tree. Each line in the file represents an edge, and the edges should be in the format: `vertex1 vertex2`. For example:

        1   2
        2   3
        2   4 
        1   3

2. **Location/Leaf labeling file** : The location/leaf labeling file contains the labels assigned to each vertex/leaf. Each line in the file corresponds to a vertex/leaf of the tree represented in the tree file, and the labels should be in the format: `vertex label`. For example:

        1   A
        2   B
        3   C

`PCCH` outputs a location labeling file, which follows the same format as described. It also returns a DOT file containing the tree along with location labeling. `PCCH` also prints `<primary location> <number of migrations> <number of comigrations> Optimal <running time (in seconds)>` on console.

### Usage

`TCC.py` takes as input a tree file and a location labeling file only, and no optional arguments.

        usage: TCC.py phylogeny location_labeling

        TCC

        positional arguments:
        phylogeny             Input phylogeny
        leaf_labeling         Input location labeling

An example execution

        $ python TCC.py t.tree l.labeling
        true         0.006915225982666

`PCC` takes as input a tree file and a location labeling file along with optional arguments.

        usage: PCC.py [-h] [-p PRIMARY] [--log] [-o OUTPUT] [-t THREADS]
                    phylogeny leaf_labeling

        PCCH

        positional arguments:
        phylogeny             Input phylogeny
        location_labeling         Input location labeling

        optional arguments:
        -h, --help            show this help message and exit
        -p PRIMARY, --primary PRIMARY
                                Primary location
        --log                 Outputs Gurobi logging
        -o OUTPUT, --output OUTPUT
                                Output folder
        -t THREADS, --threads THREADS
                            Number of threads

An example execution

        $ python PCC.py t.tree l.labeling
        a-      6       Optimal         0.366915225982666

`PCCH` takes as input a tree file and a location labeling file along with optional arguments.

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
        