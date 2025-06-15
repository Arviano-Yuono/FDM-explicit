from .node import SingleNode
from typing import List, Dict
import numpy as np

class Grid:
    def __init__(self, 
                 nodes: list[SingleNode],
                 alpha: float,
                 k: float,
                 length: float):
        assert nodes[0].type == "dirichlet" or nodes[0].type == "neumann" and nodes[-1].type == "dirichlet" or nodes[-1].type == "neumann", "First and last nodes must be dirichlet or neumann"
        
        self.nodes = nodes
        self.nx = len(nodes)
        self.dx = nodes[1].x - nodes[0].x if self.nx > 1 else 0
        self.length = length
        
        self.iteration: int = 0
        self.time: float = 0

        self.alpha: float = alpha
        self.k: float = k
        
        self.T = [node.T for node in nodes]
        self.T_history: Dict[int, Dict[float,List[float]]] = {self.iteration: {self.time: [node.T for node in nodes]}}

    def __str__(self):
        return f"Grid with {len(self)} nodes, dx={self.dx}, nx={self.nx}, T_history={self.T_history}"
    
    def __repr__(self):
        return self.__str__()
    
    def __len__(self):
        return self.nx
    
    def set_T(self, T: List[float]):
        for node_idx in range(self.nx):
            self.nodes[node_idx].T = T[node_idx]
        self.T = T

    def get_T(self, iteration: int | None = None):
        if iteration is None:
            return self.T
        else:
            return np.array(list(self.T_history[iteration].values())).reshape(self.nx, 1)
    
    def update_T_history(self, T: List[float], dt: float):
        self.time += dt
        self.iteration += 1
        for node_idx in range(self.nx):
            self.nodes[node_idx].T = T[node_idx]
        self.T_history.update({self.iteration: {self.time: T}})

    def get_T_history(self):
        return self.T_history