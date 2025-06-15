import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from fdm import Explicit, HeatSource, Bar

class HeatSimulationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("1D Heat Transfer Simulation")
        
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky="nsew")
        
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=0, column=0, sticky="nsew")
        
        bar_tab = ttk.Frame(notebook)
        source_tab = ttk.Frame(notebook)
        solver_tab = ttk.Frame(notebook)
        plot_tab = ttk.Frame(notebook)
        
        notebook.add(bar_tab, text="Bar Properties")
        notebook.add(source_tab, text="Heat Source")
        notebook.add(solver_tab, text="Solver")
        notebook.add(plot_tab, text="Plotting")
        
        self.create_bar_properties_tab(bar_tab)
        
        self.create_heat_source_tab(source_tab)
        
        self.create_solver_tab(solver_tab)
        
        self.create_plotting_tab(plot_tab)
        
        ttk.Button(main_frame, text="Run Simulation", command=self.run_simulation).grid(row=1, column=0, pady=10)

    def create_bar_properties_tab(self, parent):
        ttk.Label(parent, text="Boundary Conditions", font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)
        
        ttk.Label(parent, text="Dirichlet Boundary (T1, T2):").grid(row=1, column=0, sticky=tk.W)
        self.t1_var = tk.StringVar(value="200")
        self.t2_var = tk.StringVar(value="None")
        ttk.Entry(parent, textvariable=self.t1_var, width=10).grid(row=1, column=1, sticky=tk.W)
        ttk.Entry(parent, textvariable=self.t2_var, width=10).grid(row=1, column=2, sticky=tk.W)
        
        ttk.Label(parent, text="Neumann Boundary (dT1, dT2):").grid(row=2, column=0, sticky=tk.W)
        self.dt1_var = tk.StringVar(value="None")
        self.dt2_var = tk.StringVar(value="True")
        ttk.Entry(parent, textvariable=self.dt1_var, width=10).grid(row=2, column=1, sticky=tk.W)
        ttk.Entry(parent, textvariable=self.dt2_var, width=10).grid(row=2, column=2, sticky=tk.W)
        
        ttk.Label(parent, text="Interior Temperature:").grid(row=3, column=0, sticky=tk.W)
        self.t_interior_var = tk.StringVar(value="50")
        ttk.Entry(parent, textvariable=self.t_interior_var, width=10).grid(row=3, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Length (m):").grid(row=4, column=0, sticky=tk.W)
        self.length_var = tk.StringVar(value="0.1")
        ttk.Entry(parent, textvariable=self.length_var, width=10).grid(row=4, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Number of Nodes:").grid(row=5, column=0, sticky=tk.W)
        self.num_nodes_var = tk.StringVar(value="51")
        ttk.Entry(parent, textvariable=self.num_nodes_var, width=10).grid(row=5, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Thermal Conductivity (k):").grid(row=6, column=0, sticky=tk.W)
        self.k_var = tk.StringVar(value="237")
        ttk.Entry(parent, textvariable=self.k_var, width=10).grid(row=6, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Density (ρ):").grid(row=7, column=0, sticky=tk.W)
        self.rho_var = tk.StringVar(value="2700")
        ttk.Entry(parent, textvariable=self.rho_var, width=10).grid(row=7, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Specific Heat (c):").grid(row=8, column=0, sticky=tk.W)
        self.c_var = tk.StringVar(value="897")
        ttk.Entry(parent, textvariable=self.c_var, width=10).grid(row=8, column=1, sticky=tk.W)

    def create_heat_source_tab(self, parent):
        ttk.Label(parent, text="Heat Source Properties", font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)
        
        ttk.Label(parent, text="Heat Source Value (q):").grid(row=1, column=0, sticky=tk.W)
        self.q_value_var = tk.StringVar(value="3e5")
        ttk.Entry(parent, textvariable=self.q_value_var, width=15).grid(row=1, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Source Type:").grid(row=2, column=0, sticky=tk.W)
        self.source_type_var = tk.StringVar(value="constant")
        source_types = ["constant", "sinusoidal", "triangular_symmetric", "gaussian"]
        ttk.Combobox(parent, textvariable=self.source_type_var, values=source_types, width=15).grid(row=2, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Point Index:").grid(row=3, column=0, sticky=tk.W)
        self.point_index_var = tk.StringVar(value="25")
        ttk.Entry(parent, textvariable=self.point_index_var, width=10).grid(row=3, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Peak Position:").grid(row=4, column=0, sticky=tk.W)
        self.peak_position_var = tk.StringVar(value="0.05")
        ttk.Entry(parent, textvariable=self.peak_position_var, width=10).grid(row=4, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Frequency:").grid(row=5, column=0, sticky=tk.W)
        self.frequency_var = tk.StringVar(value="4")
        ttk.Entry(parent, textvariable=self.frequency_var, width=10).grid(row=5, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Standard Deviation:").grid(row=6, column=0, sticky=tk.W)
        self.std_var = tk.StringVar(value="0.02")
        ttk.Entry(parent, textvariable=self.std_var, width=10).grid(row=6, column=1, sticky=tk.W)

    def create_solver_tab(self, parent):
        ttk.Label(parent, text="Solver Properties", font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)
        
        ttk.Label(parent, text="Max Iterations:").grid(row=1, column=0, sticky=tk.W)
        self.max_iteration_var = tk.StringVar(value="5000")
        ttk.Entry(parent, textvariable=self.max_iteration_var, width=10).grid(row=1, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Time Step (dt):").grid(row=2, column=0, sticky=tk.W)
        self.dt_var = tk.StringVar(value="0.03")
        ttk.Entry(parent, textvariable=self.dt_var, width=10).grid(row=2, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Error Tolerance:").grid(row=3, column=0, sticky=tk.W)
        self.error_var = tk.StringVar(value="1e-6")
        ttk.Entry(parent, textvariable=self.error_var, width=10).grid(row=3, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Use Analytical Solution:").grid(row=4, column=0, sticky=tk.W)
        self.is_analytical_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(parent, variable=self.is_analytical_var).grid(row=4, column=1, sticky=tk.W)
        
        tooltip_text = "Analytical solution is only available for constant heat source with Dirichlet or Neumann BC and triangular symmetric source with Dirichlet BC"
        tooltip = ttk.Label(parent, text=tooltip_text, wraplength=200, foreground="gray")
        tooltip.grid(row=5, column=0, columnspan=2, sticky=tk.W, pady=(0, 5))

    def create_plotting_tab(self, parent):
        ttk.Label(parent, text="Plotting Properties", font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)
        
        ttk.Label(parent, text="Plot Speed:").grid(row=1, column=0, sticky=tk.W)
        self.plot_speed_var = tk.StringVar(value="25")
        ttk.Entry(parent, textvariable=self.plot_speed_var, width=10).grid(row=1, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Bar Height:").grid(row=2, column=0, sticky=tk.W)
        self.bar_height_var = tk.StringVar(value="50")
        ttk.Entry(parent, textvariable=self.bar_height_var, width=10).grid(row=2, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Number of Text Labels:").grid(row=3, column=0, sticky=tk.W)
        self.num_text_var = tk.StringVar(value="5")
        ttk.Entry(parent, textvariable=self.num_text_var, width=10).grid(row=3, column=1, sticky=tk.W)
        
        ttk.Label(parent, text="Save Plot:").grid(row=4, column=0, sticky=tk.W)
        self.save_plot_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(parent, variable=self.save_plot_var).grid(row=4, column=1, sticky=tk.W)

    def parse_boundary_conditions(self):
        #Dirichlet boundary
        t1 = float(self.t1_var.get()) if self.t1_var.get().lower() != "none" else None
        t2 = float(self.t2_var.get()) if self.t2_var.get().lower() != "none" else None
        
        #Neumann boundary
        dt1 = True if self.dt1_var.get().lower() == "true" else None
        dt2 = True if self.dt2_var.get().lower() == "true" else None
        
        return [t1, t2], [dt1, dt2]

    def run_simulation(self):
        try:
            dirichlet_boundary, neumann_boundary = self.parse_boundary_conditions()
            
            params = {
                'Dirichlet_boundary': dirichlet_boundary,
                'Neumann_boundary': neumann_boundary,
                'T_interior': float(self.t_interior_var.get()),
                'length': float(self.length_var.get()),
                'num_nodes': int(self.num_nodes_var.get()),
                'k': float(self.k_var.get()),
                'rho': float(self.rho_var.get()),
                'c': float(self.c_var.get()),
                'q_value': float(self.q_value_var.get()),
                'source_type': self.source_type_var.get(),
                'point_index': float(self.point_index_var.get()),
                'peak_position': float(self.peak_position_var.get()),
                'frequency': float(self.frequency_var.get()),
                'std': float(self.std_var.get()),
                'max_iteration': int(self.max_iteration_var.get()),
                'dt': float(self.dt_var.get()),
                'error': float(self.error_var.get()),
                'is_analytical': self.is_analytical_var.get(),
                'plot_speed': int(self.plot_speed_var.get()),
                'bar_height': int(self.bar_height_var.get()),
                'num_text': int(self.num_text_var.get()),
                'save_plot': self.save_plot_var.get()
            }
            
            alpha = np.sqrt(params['rho'] * params['c'] / params['k'])
            
            bar = Bar(params['length'], 
                     params['num_nodes'],
                     params['Dirichlet_boundary'],
                     params['Neumann_boundary'],
                     params['k'],
                     params['T_interior'],
                     alpha)
            
            analytical_grid = bar.make_grid()
            numerical_grid = bar.make_grid()
            
            heat_source = HeatSource(numerical_grid,
                                   source_type=params['source_type'],
                                   source_value=params['q_value'],
                                   source_position=params['peak_position'],
                                   std_dev=params['std'],
                                   point_index=int(params['point_index']),
                                   frequency=params['frequency'])
            heat_source.plot()

            solver = Explicit(numerical_grid,
                            params['dt'],
                            params['error'],
                            params['max_iteration'],
                            heat_source)
            
            solver.solve()
            bar.set_numerical_grid(numerical_grid)
            
            if params['is_analytical']:
                from fdm import AnalyticalSolution
                analytical_solution = AnalyticalSolution(grid=analytical_grid,
                                                      bar=bar,
                                                      T_interior=params['T_interior'],
                                                      heat_source=heat_source)
                analytical_solution.solve(verbose=True)
                bar.set_analytical_grid(analytical_solution.grid)
            
            messagebox.showinfo("Success", "Simulation completed successfully!")
            bar.plot_1d_animation(plot_speed = params['plot_speed'], analytical_solution = params['is_analytical'])
            bar.plot_2d_animation(plot_speed = params['plot_speed'], height = params['bar_height'], num_text = params['num_text'])

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = HeatSimulationGUI(root)
    root.mainloop() 