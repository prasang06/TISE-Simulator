import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

hbar_default = 1.0

# Helper: mapping of safe names for custom potential evaluation
SAFE_GLOBALS = {'__builtins__': {}}
SAFE_LOCAL_NAMES = {
    'np': np,
    'sin': np.sin,
    'cos': np.cos,
    'tan': np.tan,
    'arcsin': np.arcsin,
    'arccos': np.arccos,
    'arctan': np.arctan,
    'sinh': np.sinh,
    'cosh': np.cosh,
    'tanh': np.tanh,
    'exp': np.exp,
    'log': np.log,
    'log10': np.log10,
    'sqrt': np.sqrt,
    'abs': np.abs,
    'pi': np.pi,
    'e': np.e,
    'where': np.where
}

class TISEApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('TISE Visualizer')
        self.geometry('1200x720')
        self.style = ttk.Style(self)
        try:
            self.style.theme_use('clam')
        except Exception:
            pass

        self._build_controls()
        self._build_plot_area()

        self.x = None
        self.energies = None
        self.eigvecs = None

    def _build_controls(self):
        ctrl = ttk.Frame(self)
        ctrl.pack(side='left', fill='y', padx=10, pady=10)

        ttk.Label(ctrl, text='Potential:').grid(row=0, column=0, sticky='w')
        self.pot_var = tk.StringVar(value='Infinite Well')
        pot_cb = ttk.Combobox(
            ctrl,
            textvariable=self.pot_var,
            values=['Infinite Well', 'Finite Well', 'Harmonic Oscillator', 'Custom'],
            state='readonly'
        )
        pot_cb.grid(row=0, column=1, sticky='ew')
        pot_cb.bind('<<ComboboxSelected>>', lambda e: self._update_potential_inputs())

        ttk.Separator(ctrl, orient='horizontal').grid(
            row=1, column=0, columnspan=2, sticky='ew', pady=6)

        ttk.Label(ctrl, text='Mass (m)').grid(row=2, column=0, sticky='w')
        self.m_entry = ttk.Entry(ctrl); self.m_entry.grid(row=2, column=1); self.m_entry.insert(0, '1.0')

        ttk.Label(ctrl, text='ħ (hbar)').grid(row=3, column=0, sticky='w')
        self.hbar_entry = ttk.Entry(ctrl); self.hbar_entry.grid(row=3, column=1); self.hbar_entry.insert(0, str(hbar_default))

        ttk.Label(ctrl, text='Length (L)').grid(row=4, column=0, sticky='w')
        self.L_entry = ttk.Entry(ctrl); self.L_entry.grid(row=4, column=1); self.L_entry.insert(0, '10.0')

        ttk.Label(ctrl, text='Grid points (N)').grid(row=5, column=0, sticky='w')
        self.N_entry = ttk.Entry(ctrl); self.N_entry.grid(row=5, column=1); self.N_entry.insert(0, '800')

        self.pot_frame = ttk.Frame(ctrl)
        self.pot_frame.grid(row=6, column=0, columnspan=2, pady=(8,0), sticky='ew')
        self._update_potential_inputs()

        ttk.Separator(ctrl, orient='horizontal').grid(row=7, column=0, columnspan=2, sticky='ew', pady=6)

        ttk.Label(ctrl, text='# eigenstates to compute').grid(row=8, column=0, sticky='w')
        self.num_states_spin = ttk.Spinbox(ctrl, from_=1, to=200, width=5); self.num_states_spin.grid(row=8, column=1); self.num_states_spin.set('6')

        compute_btn = ttk.Button(ctrl, text='Compute', command=self.compute)
        compute_btn.grid(row=9, column=0, columnspan=2, pady=6, sticky='ew')

        ttk.Label(ctrl, text='Select state to display').grid(row=10, column=0, sticky='w')
        self.state_spin = ttk.Spinbox(ctrl, from_=1, to=200, width=5); self.state_spin.grid(row=10, column=1); self.state_spin.set('1')

        show_btn = ttk.Button(ctrl, text='Show state', command=self.show_state)
        show_btn.grid(row=11, column=0, columnspan=2, sticky='ew', pady=4)

        ttk.Label(ctrl, text='Energy eigenvalues').grid(row=12, column=0, columnspan=2, sticky='w', pady=(8,0))
        self.energy_list = tk.Text(ctrl, width=40, height=12)
        self.energy_list.grid(row=13, column=0, columnspan=2, pady=4)

        save_btn = ttk.Button(ctrl, text='Save plot as PNG', command=self.save_plot)
        save_btn.grid(row=14, column=0, columnspan=2, sticky='ew', pady=(10,0))

        ttk.Label(ctrl, text='Custom potential help: use Python expression in x, e.g. sin(2*pi*x/L) or np.where(x<L/2, -50, 0)').grid(row=15, column=0, columnspan=2, sticky='w')

        for child in ctrl.winfo_children():
            child.grid_configure(padx=4, pady=4)

    def _clear_pot_frame(self):
        for c in self.pot_frame.winfo_children():
            c.destroy()

    def _update_potential_inputs(self):
        self._clear_pot_frame()
        pot = self.pot_var.get()
        if pot == 'Infinite Well':
            ttk.Label(self.pot_frame, text='Infinite well: V=0 in [0,L]').grid(row=0, column=0, sticky='w')
        elif pot == 'Finite Well':
            ttk.Label(self.pot_frame, text='Well depth (positive)').grid(row=0, column=0, sticky='w')
            self.depth_entry = ttk.Entry(self.pot_frame); self.depth_entry.grid(row=0, column=1); self.depth_entry.insert(0, '50.0')
        elif pot == 'Harmonic Oscillator':
            ttk.Label(self.pot_frame, text='Omega (ω)').grid(row=0, column=0, sticky='w')
            self.omega_entry = ttk.Entry(self.pot_frame); self.omega_entry.grid(row=0, column=1); self.omega_entry.insert(0, '1.0')
        elif pot == 'Custom':
            ttk.Label(self.pot_frame, text='Custom potential (python expression in x)').grid(row=0, column=0, sticky='w')
            self.custom_entry = ttk.Entry(self.pot_frame, width=40)
            self.custom_entry.grid(row=0, column=1, sticky='ew')
            self.custom_entry.insert(0, '0.5*(x - L/2)**2')

    def _build_plot_area(self):
        plot_frame = ttk.Frame(self)
        plot_frame.pack(side='right', fill='both', expand=True)

        self.fig = Figure(figsize=(8,6), dpi=110)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('|ψ|^2 / V')
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)

        bottom = ttk.Frame(plot_frame)
        bottom.pack(fill='x')
        ttk.Label(bottom, text='Status:').pack(side='left')
        self.status_label = ttk.Label(bottom, text='Ready')
        self.status_label.pack(side='left', padx=6)

    def compute(self):
        try:
            m = float(self.m_entry.get())
            hbar = float(self.hbar_entry.get())
            L = float(self.L_entry.get())
            N = int(self.N_entry.get())
            num_states = int(self.num_states_spin.get())
            if N < 10:
                raise ValueError('N too small')
        except Exception:
            messagebox.showerror('Input error', 'Please check numeric inputs')
            return

        self.status_label.config(text='Building Hamiltonian...')
        self.update_idletasks()

        x = np.linspace(0, L, N)
        dx = x[1] - x[0]

        main_diag = -2.0 * np.ones(N)
        off_diag = 1.0 * np.ones(N-1)
        lap = (np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)) / (dx**2)
        T = - (hbar**2) / (2.0 * m) * lap

        pot = np.zeros_like(x)
        pot_type = self.pot_var.get()
        if pot_type == 'Infinite Well':
            pot[:] = 0.0
            pot[0] = pot[-1] = 1e12
        elif pot_type == 'Finite Well':
            depth = float(self.depth_entry.get())
            half = L/2
            width = L/2
            pot[:] = 0.0
            pot[(x >= (half - width/2)) & (x <= (half + width/2))] = -abs(depth)
        elif pot_type == 'Harmonic Oscillator':
            omega = float(self.omega_entry.get())
            center = L/2
            pot = 0.5 * m * (omega**2) * (x - center)**2
        elif pot_type == 'Custom':
            expr = (self.custom_entry.get() or '').strip()
            safe_dict = dict(SAFE_LOCAL_NAMES)
            safe_dict.update({'x': x, 'm': m, 'L': L, 'hbar': hbar})
            try:
                pot_val = eval(expr, SAFE_GLOBALS, safe_dict)
            except Exception as e:
                messagebox.showerror('Custom potential error', f'Error evaluating expression: {e}')
                return
            if np.isscalar(pot_val):
                pot[:] = float(pot_val)
            else:
                pot = np.array(pot_val, dtype=float)
                if pot.shape != x.shape:
                    messagebox.showerror('Custom potential error', 'Expression must return scalar or array of same length as x')
                    return

        H = T + np.diag(pot)

        try:
            eigvals, eigvecs = np.linalg.eigh(H)
        except Exception as e:
            messagebox.showerror('Solver error', str(e))
            return

        self.x = x
        self.energies = eigvals[:num_states]
        self.eigvecs = eigvecs[:, :num_states]

        self.energy_list.delete('1.0', 'end')
        for i, E in enumerate(self.energies, start=1):
            self.energy_list.insert('end', f'{i}: {E:.6g}\n')

        self.status_label.config(text='Computed')
        self.show_state()

    def show_state(self):
        if self.eigvecs is None:
            messagebox.showinfo('No data', 'Compute eigenstates first')
            return
        idx = int(self.state_spin.get()) - 1
        if idx < 0 or idx >= self.eigvecs.shape[1]:
            messagebox.showerror('Index error', 'State index out of range')
            return

        psi = self.eigvecs[:, idx]
        prob = np.real(psi * np.conj(psi))

        self.ax.clear()
        self.ax.plot(self.x, prob, label=f'|ψ_{idx+1}|^2')

        scaled_pot = self._get_current_potential()
        if scaled_pot is not None:
            pv = scaled_pot - scaled_pot.min()
            if pv.max() != 0:
                pv = pv / pv.max() * np.max(prob) * 0.9
                self.ax.plot(self.x, pv, linestyle='--', label='scaled V(x)')

        self.ax.set_xlabel('x')
        self.ax.set_ylabel('Probability density')
        self.ax.set_title(f'State {idx+1}, Energy = {self.energies[idx]:.6g}')
        self.ax.legend()
        self.ax.grid(True)
        self.canvas.draw_idle()

    def _get_current_potential(self):
        if self.x is None:
            return None
        pot_type = self.pot_var.get()
        x = self.x
        L = float(self.L_entry.get())
        if pot_type == 'Infinite Well':
            v = np.zeros_like(x)
            v[0] = v[-1] = 1.0
            return v
        elif pot_type == 'Finite Well':
            depth = float(self.depth_entry.get())
            half = L/2
            width = L/2
            v = np.zeros_like(x)
            v[(x >= (half - width/2)) & (x <= (half + width/2))] = -abs(depth)
            return v
        elif pot_type == 'Harmonic Oscillator':
            omega = float(self.omega_entry.get())
            m = float(self.m_entry.get())
            center = L/2
            return 0.5 * m * (omega ** 2) * (x - center) ** 2
        elif pot_type == 'Custom':
            try:
                expr = (self.custom_entry.get() or '').strip()
            except Exception:
                return None
            safe_dict = dict(SAFE_LOCAL_NAMES)
            safe_dict.update({'x': x, 'm': float(self.m_entry.get()), 'L': L, 'hbar': float(self.hbar_entry.get())})
            try:
                pot_val = eval(expr, SAFE_GLOBALS, safe_dict)
            except Exception:
                return None
            if np.isscalar(pot_val):
                return np.zeros_like(x) + float(pot_val)
            else:
                pv = np.array(pot_val, dtype=float)
                if pv.shape == x.shape:
                    return pv
                return None
        return None

    def save_plot(self):
        try:
            self.fig.savefig('tise_plot.png')
            messagebox.showinfo('Saved', 'Saved as tise_plot.png')
        except Exception as e:
            messagebox.showerror('Save error', str(e))

if __name__ == '__main__':
    app = TISEApp()
    app.mainloop()
