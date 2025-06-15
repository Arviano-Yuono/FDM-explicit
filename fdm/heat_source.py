from matplotlib import pyplot as plt
import numpy as np
from fdm.grid import Grid

class HeatSource:
    def __init__(self, 
                 grid: Grid, 
                 source_type: str | None = None, 
                 source_value: float | None = None, 
                 source_position: float | None = None, 
                 std_dev: float | None = None,
                 point_index: int | None = None,
                 frequency: float | None = None):
        
        self.grid = grid
        self.source_value = 0.0
        self.source_position = source_position
        self.std_dev = std_dev
        self.frequency = frequency
        self.point_index = point_index
        self.fig = None
        self.ax = None
        if source_type is not None:
            self.source_title = " ".join([text[0].upper() + text[1:] for text in list(source_type.split("_"))])
            self.source_type = source_type.lower()
        else:
            self.source_title = None
            self.source_type = None
        
        if source_type is not None:
            assert source_value is not None, "Source value is required for source type"
            self.source_value = float(source_value)

    def get_source_value(self, node_idx: int):
        if self.source_type == "point":
            assert self.point_index is not None, "Point index is required for point source"
            return self.source_value * (node_idx == self.point_index)
        elif self.source_type == "constant":
            return self.source_value
        elif self.source_type == "line":
            assert self.source_position is not None, "Source position is required for line source"
            return self.source_value * (self.source_position - self.grid.nodes[node_idx].x)
        elif self.source_type == "gaussian":
            assert self.source_position is not None, "Source position is required for gaussian source"
            assert self.std_dev is not None, "Standard deviation is required for gaussian source"
            return self.source_value * np.exp(-(self.grid.nodes[node_idx].x - self.source_position) ** 2 / (2 * self.std_dev ** 2))
        elif self.source_type == "sinusoidal":
            assert self.frequency is not None, "Frequency is required for sinusoidal source"
            return self.source_value/2 * np.sin(self.frequency * np.pi * self.grid.nodes[node_idx].x / self.grid.length) + self.source_value/2
        elif self.source_type == "triangular_symmetric":
            x = self.grid.nodes[node_idx].x
            L = self.grid.length
            peak = L/2
            return (2 * self.source_value / L) * \
                (x if x < peak else L - x)
        elif self.source_type == "square":
            assert self.source_position is not None, "Source position is required for square source"
            return self.source_value * (self.grid.nodes[node_idx].x - self.source_position)
        elif self.source_type == None:
            return 0.0
        else:
            raise ValueError(f"Invalid source type: {self.source_type}")
        
    def get_source_value_array(self):
        return np.array([self.get_source_value(node_idx) for node_idx in range(len(self.grid))])
    
    def plot(self):
        x = [node.x for node in self.grid.nodes]
        
        # Create figure and axis if they don't exist
        if self.fig is None:
            self.fig, self.ax = plt.subplots(figsize=(10, 5))
        
        # Clear the axis
        self.ax.clear()
        
        # Plot the data
        self.ax.plot(x, self.get_source_value_array())
        self.ax.set_xlabel(r"Length (m)")
        self.ax.set_ylabel(r"Heat Generation Rate ($\mathrm{W/m^3}$)")
        self.ax.set_title(rf"{self.source_title} Heat Source")
        self.ax.set_ylim(0, 1.1*self.source_value)
        self.ax.set_xlim(0, self.grid.length)
        self.ax.grid(True)
        
        plt.draw()
        plt.pause(0.001)