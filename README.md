# 1D Heat Transfer Simulation

A Python-based simulation tool for analyzing one-dimensional heat transfer problems using the Finite Difference Method (FDM).

## Features

- 1D heat transfer simulation using explicit FDM
- Support for both Dirichlet and Neumann boundary conditions
- Configurable heat source, with limited support of analytical comparison of:
  - constant heat source with 2 Dirchlet BC
  - constant heat source with Dirichlet and Neumann BC
  - triangle heat source with 2 Dirichlet BC
- Interactive GUI for parameter input
- Real-time visualization of temperature distribution

## Requirements

- Python 3.x
- numpy
- matplotlib
- tkinter
- tqdm

## Installation

1. Clone the repository:
   `git pull https://github.com/Arviano-Yuono/FDM-explicit.git`
2. Go to the repo dir
   `cd FDM-explicit`
3. Install the dependencies
   `pip install -r requirements.txt`
4. Run `main.py` using
   `python main.py`
