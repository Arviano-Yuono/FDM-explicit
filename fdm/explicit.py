from .grid import Grid
from tqdm import tqdm
from .heat_source import HeatSource
class Explicit:
    def __init__(self, 
                 grid: Grid, 
                 dt: float, 
                 error: float, 
                 max_iteration: int,
                 heat_source: HeatSource):
        
        self.grid = grid
        self.dt = dt        
        self.error = error
        self.max_iteration = max_iteration
        self.heat_source = heat_source

        # grid properties
        self.dx = self.grid.dx
        self.alpha = self.grid.alpha
        self.k = self.grid.k


    def solve(self, verbose: bool = True):
        if verbose:
            print("Solving...")
        loop = tqdm(range(1, self.max_iteration), disable=not verbose)
        for iteration in loop:
            self.sum_error = 0
            for node_idx in range(len(self.grid)):
                # print(self.grid.nodes[node_idx])
                if self.grid.nodes[node_idx].type == "interior":
                    self.grid.nodes[node_idx].T = self.get_next_T(node_idx)
                    self.sum_error += self.get_error(node_idx)
                elif self.grid.nodes[node_idx].type == "neumann":
                    closest_node_idx = node_idx + 1 if node_idx==0 else node_idx - 1
                    self.grid.nodes[node_idx].T = self.get_next_T(closest_node_idx)
                    self.sum_error += self.get_error(closest_node_idx)
            if self.sum_error / len(self.grid) < self.error:
                print(f"Converged at iteration {iteration}")
                break
            self.grid.update_T_history(T = [node.T for node in self.grid.nodes], dt = self.dt)
            loop.set_postfix(error=self.sum_error / len(self.grid))

    def get_next_T(self, node_idx: int):
        Ti = self.grid.nodes[node_idx].T
        Ti_next = self.grid.nodes[node_idx+1].T
        Ti_prev = self.grid.nodes[node_idx-1].T
        C = self.dt/((self.alpha**2)*(self.dx**2))
        heat_source_term = self.heat_source.get_source_value(node_idx) * self.dt / (self.k * self.alpha**2)
        Ti_new = Ti + C * (Ti_next - 2 * Ti + Ti_prev) + heat_source_term
        return Ti_new
        
    def get_error(self, node_idx: int):
        return abs(self.grid.nodes[node_idx].T - self.get_next_T(node_idx))
    
    def get_max_error(self):
        return max(self.get_error(node_idx) for node_idx in range(len(self.grid)))
    
    def is_converged(self):
        return self.sum_error / len(self.grid) < self.error