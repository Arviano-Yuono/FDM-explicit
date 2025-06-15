from .grid import Grid
from .node import SingleNode
from IPython.display import clear_output
import matplotlib.pyplot as plt
import numpy as np
from typing import Union, List

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
        self.fig_1d = None
        self.fig_2d = None
        self.ax_1d = None
        self.ax_2d = None
        self.im_2d = None
        
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

    def plot_2d(self, plot_speed: int = 10, height: int = 10, num_text: int = 5):
        assert self.numerical_grid is not None, "Numerical grid is not set"

        # Calculate text positions to be evenly spread
        text_positions = np.linspace(0, self.num_nodes - 1, num_text + 2, dtype=int)[1:-1]
        
        # Create figure and axis if they don't exist
        if self.fig_2d is None:
            self.fig_2d, self.ax_2d = plt.subplots(figsize=(12, 2))
            self.cbar = None
        
        # Get all iterations including the last one
        iterations = list(range(0, len(self.numerical_grid.get_T_history())-1, plot_speed))
        iterations.append(len(self.numerical_grid.get_T_history())-1)
            
        for iteration in iterations:
            clear_output(wait=True)
            time, last_temp = list(self.numerical_grid.get_T_history()[iteration].items())[0]
            
            temp_1d = np.array(last_temp)
            temp_2d = np.tile(temp_1d, (height, 1))  # fixed small height for bar look
            vmin = np.min(temp_1d)
            vmax = np.max(temp_1d)
            
            # Clear the axis
            self.ax_2d.clear()
            
            # Create or update the image
            self.im_2d = self.ax_2d.imshow(temp_2d, 
                                          aspect='auto',
                                          cmap='hot',
                                          extent=(0.0, float(self.length), 0.0, 1.0),  # y from 0 to 1
                                          interpolation='nearest',
                                          vmin=vmin,
                                          vmax=vmax)
            # Update or create colorbar
            if self.cbar is not None:
                try:
                    self.cbar.remove()
                except Exception:
                    pass
                self.cbar = None
            self.cbar = self.fig_2d.colorbar(self.im_2d, ax=self.ax_2d, label=r'$\mathrm{Temperature\ (K)}$')
            
            # Set labels and title
            self.ax_2d.set_xlabel(r'$\mathrm{Position\ (m)}$')
            self.ax_2d.set_ylabel(r'$\mathrm{Normalized\ Height}$')
            self.ax_2d.set_title(r'$\mathrm{Temperature\ Distribution\ -\ Iteration\ ' + str(iteration) + ' (t = ' + f"{time:.4f}" + ' s)}$')
            self.ax_2d.set_ylim(0, 1)
            self.ax_2d.set_xlim(0, float(self.length))
            self.ax_2d.set_yticks([])  # Optional: remove y-ticks for clarity
            
            # Add temperature labels at evenly spaced positions
            for j in text_positions:
                temp = temp_1d[j]
                self.ax_2d.text(j * self.dx, 0.5, f'{float(temp):.1f}°C', 
                              ha='center', va='center', color='green', fontsize=12)
            
            plt.draw()
            plt.pause(0.001)

    def plot_1d_animation(self, plot_speed: int = 10, analytical_solution = False):
        assert self.numerical_grid is not None, "Numerical grid is not set"
        self.T_boundary = [self.numerical_grid.nodes[0].T, self.numerical_grid.nodes[-1].T]
        x = np.array([node.x for node in self.numerical_grid.nodes])
        
        # Create figure and axis if they don't exist
        if self.fig_1d is None:
            self.fig_1d, self.ax_1d = plt.subplots(figsize=(10, 5))
        
        # Get all iterations including the last one
        iterations = list(range(0, len(self.numerical_grid.get_T_history())-1, plot_speed))
        iterations.append(len(self.numerical_grid.get_T_history())-1)
            
        for iteration in iterations:
            clear_output(wait=True)
            time, numerical_last_temp = list(self.numerical_grid.get_T_history()[iteration].items())[0]
            
            # Clear the axis
            self.ax_1d.clear()
            
            # Set labels and title
            self.ax_1d.set_xlabel(r"$\mathrm{Length\ (m)}$")
            self.ax_1d.set_ylabel(r"$\mathrm{Temperature\ (K)}$")
            self.ax_1d.set_title(r"$\mathrm{Temperature\ Distribution\ -\ Iteration\ " + str(iteration) + " (t = " + f"{time:.4f}" + " s)}$")
            self.ax_1d.set_xlim(0, max(x))
            self.ax_1d.grid(True, linestyle='--', alpha=0.5)
            self.ax_1d.set_ylim(0, 1.1*max(numerical_last_temp))
            
            if analytical_solution:
                analytical_last_temp = self.analytical_grid.get_T()
                self.ax_1d.plot(x, analytical_last_temp, label=r'$\mathrm{Analytical\ Solution}$', linewidth=2.5)
                self.ax_1d.plot(x, numerical_last_temp, label=r'$\mathrm{Numerical\ Solution}$', linestyle='--', color='red', linewidth=2)
                self.ax_1d.set_ylim(0, 1.1*max(max(analytical_last_temp), max(numerical_last_temp)))
            else:
                self.ax_1d.plot(x, numerical_last_temp, label=r'$\mathrm{Numerical\ Solution}$', color='red', linewidth=2)

            self.ax_1d.legend()
            plt.draw()
            plt.pause(0.0001)

    def plot_2d_animation(self, plot_speed: int = 10, height: int = 10, num_text: int = 5):
        import matplotlib
        assert self.numerical_grid is not None, "Numerical grid is not set"

        def setup_figure():
            if not hasattr(self, 'fig_2d_anim') or self.fig_2d_anim is None:
                self.fig_2d_anim, self.ax_2d_anim = plt.subplots(figsize=(12, 2), constrained_layout=True)
                self.im_2d_anim = None
                self.cbar_anim = None

        def draw_frame(temp_1d, iteration, time, first_frame=False):
            temp_2d = np.tile(temp_1d, (height, 1))
            vmin, vmax = np.min(temp_1d), np.max(temp_1d)
            text_positions = np.linspace(0, self.num_nodes - 1, num_text + 2, dtype=int)[1:-1]

            if first_frame or self.im_2d_anim is None:
                self.ax_2d_anim.clear()
                self.im_2d_anim = self.ax_2d_anim.imshow(
                    temp_2d, aspect='auto', cmap='hot',
                    extent=(0.0, float(self.length), 0.0, 1.0),
                    interpolation='nearest', vmin=vmin, vmax=vmax
                )
                if self.cbar_anim is not None:
                    try: self.cbar_anim.remove()
                    except Exception: pass
                self.cbar_anim = self.fig_2d_anim.colorbar(self.im_2d_anim, ax=self.ax_2d_anim, label=r'$\mathrm{Temperature\ (K)}$', pad=0.02)
                self.ax_2d_anim.set(
                    xlabel=r'$\mathrm{Position\ (m)}$',
                    ylabel=r'$\mathrm{Normalized\ Height}$',
                    title=fr'$\mathrm{{Temperature\ Distribution\ -\ Iteration\ {iteration} (t = {time:.4f} s)}}$',
                    ylim=(0, 1), xlim=(0, float(self.length))
                )
                self.ax_2d_anim.set_yticks([])
            else:
                self.im_2d_anim.set_data(temp_2d)
                self.im_2d_anim.set_clim(vmin, vmax)
                self.cbar_anim.update_normal(self.im_2d_anim)
                self.ax_2d_anim.set_title(fr'$\mathrm{{Temperature\ Distribution\ -\ Iteration\ {iteration} (t = {time:.4f} s)}}$')

            # Remove old text labels if any
            if hasattr(self, 'text_labels_anim'):
                for txt in self.text_labels_anim:
                    txt.remove()
            self.text_labels_anim = []
            for j in text_positions:
                txt = self.ax_2d_anim.text(j * self.dx, 0.5, f'{float(temp_1d[j]):.1f}°C',
                                           ha='center', va='center', color='green', fontsize=12)
                self.text_labels_anim.append(txt)
            plt.draw()
            plt.pause(0.001)

        setup_figure()
        iterations = list(range(0, len(self.numerical_grid.get_T_history())-1, plot_speed))
        iterations.append(len(self.numerical_grid.get_T_history())-1)
        first = True
        for iteration in iterations:
            clear_output(wait=True)
            time, last_temp = list(self.numerical_grid.get_T_history()[iteration].items())[0]
            temp_1d = np.array(last_temp)
            draw_frame(temp_1d, iteration, time, first_frame=first)
            first = False