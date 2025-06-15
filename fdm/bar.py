from .grid import Grid
from .node import SingleNode
from IPython.display import clear_output
import matplotlib.pyplot as plt
import numpy as np
from typing import Union, List
import time

class Bar:
    def __init__(self, 
                 length: float, 
                 num_nodes: int,
                 dirichlet_boundary: List[Union[float, None]],
                 neumann_boundary: List[Union[float, None]],
                 k: float, 
                 T_interior: float, 
                 alpha: float):
        self.length = length
        self.k = k
        self.T_interior = T_interior
        self.dirichlet_boundary = dirichlet_boundary
        self.neumann_boundary = neumann_boundary
        self.num_nodes = num_nodes
        self.dx = self.length / (self.num_nodes - 1)
        self.alpha = alpha
        
    def make_grid(self):
        nodes = []
        # first node
        if self.dirichlet_boundary[0] is not None:
            nodes.append(SingleNode(name = f"node_0", x = 0, T = self.dirichlet_boundary[0], type = "dirichlet"))
        elif self.neumann_boundary[0] is not None:
            nodes.append(SingleNode(name = f"node_0", x = 0, T = self.T_interior, type = "neumann"))
        elif self.dirichlet_boundary[0] is not None and self.neumann_boundary[0] is not None:
            nodes.append(SingleNode(name = f"node_0", x = 0, T = self.dirichlet_boundary[0], type = "dirichlet_neumann"))
        else:
            raise ValueError("Dirichlet or Neumann boundary is required")
        
        # interior nodes
        for i in range(1, self.num_nodes - 1):
            nodes.append(SingleNode(name = f"node_{i}", x = i * self.dx, T = self.T_interior, type = "interior"))
        
        # last node
        if self.dirichlet_boundary[1] is not None:
            nodes.append(SingleNode(name = f"node_{self.num_nodes - 1}", x = self.length, T = self.dirichlet_boundary[1], type = "dirichlet"))
        elif self.neumann_boundary[1] is not None:
            nodes.append(SingleNode(name = f"node_{self.num_nodes - 1}", x = self.length, T = self.T_interior, type = "neumann"))
        elif self.dirichlet_boundary[1] is not None and self.neumann_boundary[1] is not None:
            nodes.append(SingleNode(name = f"node_{self.num_nodes - 1}", x = self.length, T = self.dirichlet_boundary[1], type = "dirichlet_neumann"))
        else:
            raise ValueError("Dirichlet or Neumann boundary is required")
        
        return Grid(nodes, 
                    alpha = self.alpha, 
                    k = self.k,
                    length=self.length)
    
    def set_numerical_grid(self, numerical_grid):
        self.numerical_grid = numerical_grid
    
    def set_analytical_grid(self, analytical_grid):
        self.analytical_grid = analytical_grid

    def plot_1d(self, plot_speed: int = 10, analytical_solution = False):
        assert self.numerical_grid is not None, "Numerical grid is not set"
        self.T_boundary = [self.numerical_grid.nodes[0].T, self.numerical_grid.nodes[-1].T]
        x = np.array([node.x for node in self.numerical_grid.nodes])
        
        plt.figure(figsize=(10, 5))
        # Get all iterations including the last one
        iterations = list(range(0, len(self.numerical_grid.get_T_history())-1, plot_speed))
        iterations.append(len(self.numerical_grid.get_T_history())-1)
            
        for iteration in iterations:
            clear_output(wait=True)
            time, numerical_last_temp = list(self.numerical_grid.get_T_history()[iteration].items())[0]
            plt.clf()
            plt.xlabel(r"$\mathrm{Length\ (m)}$")
            plt.ylabel(r"$\mathrm{Temperature\ (K)}$")
            plt.title(r"$\mathrm{Temperature\ Distribution\ -\ Iteration\ " + str(iteration) + " (t = " + f"{time:.4f}" + " s)}$")
            plt.xlim(0, max(x))
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.ylim(0, 1.1*max(numerical_last_temp))
            
            if analytical_solution:
                analytical_last_temp = self.analytical_grid.get_T()
                plt.plot(x, analytical_last_temp, label=r'$\mathrm{Analytical\ Solution}$', linewidth=2.5)
                plt.plot(x, numerical_last_temp, label=r'$\mathrm{Numerical\ Solution}$', linestyle='--', color='red', linewidth=2)
                plt.ylim(0, 1.1*max(max(analytical_last_temp), max(numerical_last_temp)))
            else:
                plt.plot(x, numerical_last_temp, label=r'$\mathrm{Numerical\ Solution}$', color='red', linewidth=2)

            plt.legend()
            plt.show()

    def plot_2d(self, plot_speed: int = 10, height: int = 50, num_text: int = 5):
        assert self.numerical_grid is not None, "Numerical grid is not set"
        self.T_boundary = [self.numerical_grid.nodes[0].T, self.numerical_grid.nodes[-1].T]
        vmin = min(self.T_boundary)
        vmax = max(self.T_boundary)

        num_text = int(self.num_nodes / num_text)
        
        # Get all iterations including the last one
        iterations = list(range(0, len(self.numerical_grid.get_T_history())-1, plot_speed))
        iterations.append(len(self.numerical_grid.get_T_history())-1)
            
        for iteration in iterations:
            clear_output(wait=True)
            time, last_temp = list(self.numerical_grid.get_T_history()[iteration].items())[0]
            
            temp_1d = np.array(last_temp)
            temp_2d = np.tile(temp_1d, (height, 1))
            
            plt.figure(figsize=(12, 4))
            plt.imshow(temp_2d, 
                      aspect='auto',
                      cmap='hot',
                      extent=(0.0, float(self.length), 0.0, 1.0),
                      interpolation='nearest',
                      vmin=vmin,
                      vmax=vmax)
            
            plt.colorbar(label=r'$\mathrm{Temperature\ (K)}$')
            plt.xlabel(r'$\mathrm{Position\ (m)}$')
            plt.ylabel(r'$\mathrm{Normalized\ Height}$')
            plt.title(r'$\mathrm{Temperature\ Distribution\ -\ Iteration\ ' + str(iteration) + ' (t = ' + f"{time:.4f}" + ' s)}$')
            
            for j, temp in enumerate(temp_1d):
                if j % num_text == 0 and j != 0 and j != temp_1d.size - 1:
                    plt.text(j * self.dx, 0.5, f'{float(temp):.1f}°C', 
                            ha='center', va='center', color='green', fontsize=12)
            
            plt.show()