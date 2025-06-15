from .grid import Grid
from .explicit import Explicit

class FDM:
    def __init__(self, grid: Grid, explicit: Explicit):
        self.grid = grid
        self.explicit = explicit

    def solve(self):
        pass