from matplotlib import pyplot as plt
from fdm import HeatSource
import numpy as np


class AnalyticalSolution:
    def __init__(self, grid, T_interior, heat_source, bar):
        self.grid = grid
        self.T_interior = T_interior
        self.heat_source = heat_source
        self.bar = bar
        
    def solve(self, verbose = False):
        T_analytical = np.zeros(self.grid.nx)
        x = np.array([node.x for node in self.grid.nodes])
        L = self.bar.length
        k = self.bar.k
        q = self.heat_source.source_value
        T1 = self.grid.nodes[0].T
        T2 = self.grid.nodes[-1].T

        if self.heat_source.source_type == "constant":
            # Check for Neumann boundary condition on the right
            if self.grid.nodes[-1].type == "neumann": 
                T_analytical = T1 - q / k * x * (x / 2 - L)
                self.grid.set_T(T_analytical)
            # Dirichlet boundary conditions on both sides
            else:
                T_analytical = -q * x**2 / (2 * k) + (T2 - T1 + q * L**2 / (2 * k)) * x / L + T1
                self.grid.set_T(T_analytical)

        elif self.heat_source.source_type == "triangular_symmetric":
            if self.grid.nodes[0].type == "dirichlet" and self.grid.nodes[-1].type == "dirichlet":
                T0 = self.grid.nodes[0].T
                T1 = self.grid.nodes[-1].T
                a = q * L**2 / k  # non-dimensional group

                # Match MATLAB's formulation
                A = np.array([
                    [0, 1, 1],
                    [-0.5, 0.5, 1],
                    [1, -1, 0]
                ])
                b = np.array([
                    T1 + 2 * a / 3,
                    T0 + a / 6,
                    0.5 * a
                ])
                C = np.linalg.solve(A, b)
                C1 = C[0] / L
                C3 = C[1] / L
                C4 = C[2]

                T_analytical = np.zeros_like(x)
                for i, xi in enumerate(x):
                    if xi <= L / 2:
                        T_analytical[i] = -q * xi**3 / (3 * k * L) + T0 + C1 * xi
                    else:
                        T_analytical[i] = q * (xi**3 / 3 - L * xi**2) / (k * L) + C3 * xi + C4

                self.grid.set_T(T_analytical)

        elif self.heat_source.source_type == "sinusoidal":
            n = self.heat_source.frequency # In the MATLAB code, n is the order, which acts as a frequency parameter
            if n is not None:
                term1 = (q * L**2 / (k * (n * np.pi)**2)) * np.sin(n * np.pi * x / L)
                term2 = (T2 - T1) / L * x
                term3 = T1
                self.T_analytical = term1 + term2 + term3
                self.grid.set_T(T_analytical)
                if verbose:
                    print(f"Sinusoidal heat source analytical solution: {self.T_analytical}")
            else:
                 raise ValueError("Frequency 'n' is required for sinusoidal heat source analytical solution.")

    def get_analytical_solution(self):
        return self.grid.get_T()

    def plot(self):
        x = np.array([node.x for node in self.grid.nodes])
        plt.figure(figsize=(10, 5))
        plt.xlabel(r"$\mathrm{Length\ (m)}$")
        plt.ylabel(r"$\mathrm{Temperature \ (K)}$")
        plt.title(r"$\mathrm{Temperature\ Distribution}$")
        plt.ylim(0, 1.1*max(self.grid.get_T()))
        plt.xlim(0, max(x))
        plt.plot(x, self.grid.get_T(), label=r'$\mathrm{Analytical\ Solution}$')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.show()
