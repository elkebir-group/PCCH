# Simulation

Phylogenetic trees with leaf labelings where MACHINA returns temporally inconsistent solutions are simulated using scripts `sim_script.py` (for disjoint cycles) and `sim_script_cmplx.py` (for complex cycles). Both the scripts are written in python and requires `Python 3.7` or newer. They also require an external python library called `networkx`, which can be installed by `pip` or `conda`.

  $ pip install networkx  
  $ conda install -c anaconda networkx  

Running the scripts returns a random instance with a phylogenetic tree and the corresponding leaf labeling. None of the scripts require any input parameter.

  $ python3 sim_script.py  
  $ python3 sim_script_cmplx.py  
