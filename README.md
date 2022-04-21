# Landscape-of-DM
## Requirements
This code is written in Python 3.7 and require <a href="https://briansimulator.org/" target="_blank">Brian2</a>.

## Reproducing results from the paper
These files reproduce all results for Fig. 2 of the paper "Leijun Ye and Chunhe Li*, Quantifying the landscape of decision making from spiking neural networks. Frontiers In Computational Neuroscience 15: 740601 (2021)."

## Usage Instructions:
1. You can start by running the Simulation.py. This will simulate the spiking neural network for decision making to get the firing rate activity for a long time. The simulations is paralleled using the multiprocessing package of Python. This file is written with the code in <a href="https://github.com/xjwanglab/book/blob/master/wang2002/wang2002.py" target="_blank">Wang Lab</a> as a reference.
2. Run the Lanscape.m to map out the potential landscape in Fig. 2(A).
3. Run the Flux.m to calculate the probablistic flux in Fig. 2(B).
4. Run the Transition_path.m to calculate the path in Fig. 2(C). The Find_region.m is used to get the spherical areas for attractors.
