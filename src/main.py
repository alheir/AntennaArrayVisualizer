import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import ttk, messagebox
from scipy.integrate import simpson  # Actualizado: simpson en lugar de simps (obsoleto)

class AntennaArrayGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Visualizador de Arreglos de Antenas Lineales")
        self.root.geometry("1400x900")
        self.root.configure(bg="#f0f0f0")
        
        # Variables por defecto
        self.N = tk.IntVar(value=4)
        self.d_lambda = tk.DoubleVar(value=0.5)
        self.ant_type = tk.StringVar(value="dipolo")
        self.plot_type = tk.StringVar(value="polar")
        
        # Amplitudes y fases por antena
        self.amplitudes = [tk.DoubleVar(value=1.0) for _ in range(10)]
        self.phases = [tk.DoubleVar(value=0.0) for _ in range(10)]
        
        # Contenedor principal horizontal
        main_container = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # ===== PANEL IZQUIERDO: Controles (25%) =====
        left_frame = ttk.Frame(main_container)
        main_container.add(left_frame, weight=1)
        
        # Título
        title_label = ttk.Label(left_frame, text="⚙ CONFIGURACIÓN", font=("Arial", 14, "bold"))
        title_label.pack(pady=10, padx=5)
        
        # Notebook para separar secciones
        notebook = ttk.Notebook(left_frame)
        notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Tab 1: Parámetros Globales
        params_frame = ttk.Frame(notebook, padding=10)
        notebook.add(params_frame, text="Parámetros")
        
        # Grid de parámetros con mejor formato
        ttk.Label(params_frame, text="Número de antenas:", font=("Arial", 10)).grid(row=0, column=0, sticky="w", pady=8)
        entry_n = ttk.Entry(params_frame, textvariable=self.N, width=15)
        entry_n.grid(row=0, column=1, sticky="ew", padx=5)
        
        ttk.Label(params_frame, text="Separación (λ):", font=("Arial", 10)).grid(row=1, column=0, sticky="w", pady=8)
        entry_d = ttk.Entry(params_frame, textvariable=self.d_lambda, width=15)
        entry_d.grid(row=1, column=1, sticky="ew", padx=5)
        
        ttk.Label(params_frame, text="Tipo de antena:", font=("Arial", 10)).grid(row=2, column=0, sticky="w", pady=8)
        combo_ant = ttk.Combobox(params_frame, textvariable=self.ant_type, 
                                  values=["isotropica", "dipolo", "monopolo"], state="readonly", width=13)
        combo_ant.grid(row=2, column=1, sticky="ew", padx=5)
        
        ttk.Label(params_frame, text="Tipo de gráfico:", font=("Arial", 10)).grid(row=3, column=0, sticky="w", pady=8)
        combo_plot = ttk.Combobox(params_frame, textvariable=self.plot_type, 
                                   values=["polar", "cartesiano"], state="readonly", width=13)
        combo_plot.grid(row=3, column=1, sticky="ew", padx=5)
        
        params_frame.columnconfigure(1, weight=1)
        
        # Botón Calcular prominente
        ttk.Button(params_frame, text="▶ CALCULAR", command=self.plot_pattern).grid(row=4, column=0, columnspan=2, pady=15, sticky="ew")
        
        # Tab 2: Configuración de Antenas
        ant_tab = ttk.Frame(notebook)
        notebook.add(ant_tab, text="Antenas")
        
        # Scrollable frame para antenas
        ant_canvas = tk.Canvas(ant_tab, bg="white", highlightthickness=0)
        scrollbar = ttk.Scrollbar(ant_tab, orient="vertical", command=ant_canvas.yview)
        scrollable_frame = ttk.Frame(ant_canvas, padding=5)
        
        scrollable_frame.bind("<Configure>", lambda e: ant_canvas.configure(scrollregion=ant_canvas.bbox("all")))
        ant_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        ant_canvas.configure(yscrollcommand=scrollbar.set)
        
        # Headers
        header_bg = "#e0e0e0"
        ttk.Label(scrollable_frame, text="Ant.", font=("Arial", 9, "bold"), background=header_bg).grid(row=0, column=0, padx=3, pady=5, sticky="ew")
        ttk.Label(scrollable_frame, text="Amplitud", font=("Arial", 9, "bold"), background=header_bg).grid(row=0, column=1, padx=3, pady=5, sticky="ew")
        ttk.Label(scrollable_frame, text="Fase (°)", font=("Arial", 9, "bold"), background=header_bg).grid(row=0, column=2, padx=3, pady=5, sticky="ew")
        
        self.ant_widgets = []
        for i in range(10):
            row = i + 1
            bg_color = "#f9f9f9" if i % 2 == 0 else "#ffffff"
            
            lbl = ttk.Label(scrollable_frame, text=f"{i+1}", background=bg_color)
            lbl.grid(row=row, column=0, padx=3, pady=3, sticky="ew")
            
            amp_entry = ttk.Entry(scrollable_frame, textvariable=self.amplitudes[i], width=8)
            amp_entry.grid(row=row, column=1, padx=3, pady=3, sticky="ew")
            
            phase_entry = ttk.Entry(scrollable_frame, textvariable=self.phases[i], width=8)
            phase_entry.grid(row=row, column=2, padx=3, pady=3, sticky="ew")
            
            self.ant_widgets.append((lbl, amp_entry, phase_entry))
        
        ant_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # ===== PANEL CENTRAL/DERECHO: Gráficos (75%) =====
        right_frame = ttk.Frame(main_container)
        main_container.add(right_frame, weight=2)
        
        # Matplotlib Figure mejorado
        self.fig = Figure(figsize=(12, 8), dpi=90)
        self.fig.patch.set_facecolor('white')
        self.canvas = FigureCanvasTkAgg(self.fig, master=right_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Cursor: Conectar evento de motion
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        
        # Frame de información en la parte inferior
        info_frame = ttk.Frame(root)
        info_frame.pack(fill="x", padx=5, pady=5)
        
        # Resultados con mejor formato
        ttk.Label(info_frame, text="📊 RESULTADOS:", font=("Arial", 10, "bold")).pack(side="left", padx=5)
        
        self.results_text = ttk.Label(info_frame, text="Directividad: -- | HPBW: -- | Cursor: --", 
                                      font=("Arial", 9), relief=tk.SUNKEN)
        self.results_text.pack(side="left", fill="x", expand=True, padx=10)
        
        self.cursor_info = ttk.Label(info_frame, text="θ: --", font=("Arial", 9), relief=tk.SUNKEN)
        self.cursor_info.pack(side="right", padx=5)
        
        # Actualizar visibilidad al cambiar N
        self.N.trace_add("write", self.update_antenna_widgets)
        
        # Plot inicial
        self.plot_pattern()
    
    def element_pattern(self, theta, ant_type):
        """Patrón elemental normalizado (|E_elem(θ)|)"""
        if ant_type == "isotropica":
            return np.ones_like(theta)
        elif ant_type == "dipolo":
            sin_theta = np.maximum(np.abs(np.sin(theta)), 1e-10)
            F = np.abs(np.cos((np.pi / 2) * np.cos(theta)) / sin_theta)
            return F / np.max(F)
        elif ant_type == "monopolo":
            F = np.where(theta <= np.pi / 2, self.element_pattern(theta, "dipolo"), 0.0)
            return F / np.max(F)
    
    def array_factor(self, theta, d_lambda, amplitudes, phases_rad):
        """Factor de arreglo para array lineal con amplitudes y fases individuales"""
        kd = 2 * np.pi * d_lambda  # k = 2π/λ, d en λ
        N = len(amplitudes)
        AF = np.zeros_like(theta, dtype=complex)
        for n in range(N):
            psi = kd * n * np.cos(theta) + phases_rad[n]
            AF += amplitudes[n] * np.exp(1j * psi)
        AF_mag = np.abs(AF)
        return AF_mag / np.max(AF_mag) if np.max(AF_mag) > 0 else np.ones_like(theta)
    
    def total_pattern(self, theta):
        """Patrón total: elemento * arreglo"""
        N = self.N.get()
        d = self.d_lambda.get()
        typ = self.ant_type.get()
        amps = [self.amplitudes[i].get() for i in range(N)]
        phases_deg = [self.phases[i].get() for i in range(N)]
        phases_rad = np.deg2rad(phases_deg)
        E_elem = self.element_pattern(theta, typ)
        AF = self.array_factor(theta, d, amps, phases_rad)
        return E_elem * AF  # Campo total normalizado
    
    def calculate_directivity(self, U_norm):
        """Directividad = 4π * U_max / ∫ U dΩ (aprox. integral numérica)"""
        # U = |E|^2, normalizado
        U = U_norm ** 2
        U_max = np.max(U)
        # Integral sobre esfera: ∫ U sinθ dθ dφ ≈ 2π ∫ U sinθ dθ (azimutal simetría)
        integ = 2 * np.pi * simpson(U * np.sin(theta), theta)  # Actualizado: simpson
        D = 4 * np.pi * U_max / integ if integ > 0 else 0
        return D
    
    def calculate_hpbw(self, U_norm):
        """Ancho de haz a media potencia (HPBW en rad, aprox. en plano principal)"""
        U = U_norm ** 2
        U_max = np.max(U)
        half = U_max / 2
        # Encontrar θ donde U(θ) = half, simétrico alrededor θ=π/2 para broadside
        idx = np.where(U > half)[0]
        if len(idx) > 0:
            theta_half = theta[idx[0]] if idx[0] > 0 else 0
            hpbw_rad = 2 * (np.pi/2 - theta_half)  # Aprox. para broadside
            return np.rad2deg(hpbw_rad)
        return 0
    
    def plot_pattern(self):
        global theta  # Para cursor
        theta = np.linspace(0, np.pi, 1000)
        U_norm = self.total_pattern(theta)
        
        # Calcular salidas
        directivity = self.calculate_directivity(U_norm)
        hpbw = self.calculate_hpbw(U_norm)
        
        # Almacenar integral para cursor
        U = U_norm ** 2
        self.integral_U = 2 * np.pi * simpson(U * np.sin(theta), theta)
        
        # Mostrar resultados en formato compacto
        results_str = f"Directividad: {directivity:.2f} | HPBW: {hpbw:.2f}°"
        self.results_text.config(text=results_str)
        
        # Limpiar figura
        self.fig.clear()
        
        # Visualización 2D
        self._plot_2d_pattern(theta, U_norm)
        
        self.fig.suptitle(f"N={self.N.get()} | d/λ={self.d_lambda.get():.2f} | {self.ant_type.get().capitalize()}", 
                         fontsize=11, fontweight='bold', y=0.995)
        self.fig.tight_layout(rect=[0, 0, 1, 0.99])
        self.canvas.draw()
    
    def _plot_2d_pattern(self, theta, U_norm):
        """Plotea el patrón en 2D: diagrama de radiación principal + disposición de antenas"""
        self.ax1 = self.fig.add_subplot(2, 1, 1, projection='polar' if self.plot_type.get() == "polar" else None)
        self.ax2 = self.fig.add_subplot(2, 1, 2)
        
        # Diagrama de radiación (plano vertical - vista más representativa)
        if self.plot_type.get() == "polar":
            self.ax1.plot(theta, U_norm, 'b-', linewidth=2.5)
            self.ax1.fill(theta, U_norm, 'b', alpha=0.2)
            self.ax1.set_title("Diagrama de Radiación (Coordenadas Polares)", fontsize=12, fontweight='bold', pad=20)
            self.ax1.grid(True, alpha=0.3)
        else:  # Cartesiano
            r = U_norm
            x = r * np.sin(theta)
            y = r * np.cos(theta)
            self.ax1.plot(x, y, 'b-', linewidth=2.5)
            self.ax1.fill(x, y, 'b', alpha=0.2)
            self.ax1.set_title("Diagrama de Radiación (Coordenadas Cartesianas)", fontsize=12, fontweight='bold', pad=10)
            self.ax1.set_aspect('equal')
            self.ax1.grid(True, alpha=0.3)
            self.ax1.set_xlabel('Eje X', fontsize=10)
            self.ax1.set_ylabel('Ángulo θ (grados)', fontsize=10)
            
            # Configurar eje Y en grados, centrado en 0°
            theta_deg = np.rad2deg(theta)
            y_ticks = np.array([0, 30, 60, 90, 120, 150, 180])
            self.ax1.set_yticks(np.cos(np.deg2rad(y_ticks)))
            self.ax1.set_yticklabels([f'{int(t)}°' for t in y_ticks])
            
            # Centrar en 0° (que corresponde a y=1)
            self.ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.3, linewidth=1)
        
        # Diagrama de disposición de antenas
        N = self.N.get()
        d = self.d_lambda.get()
        amps = [self.amplitudes[i].get() for i in range(N)]
        phases = [self.phases[i].get() for i in range(N)]
        positions = [i * d for i in range(N)]
        sizes = [max(amp * 150, 50) for amp in amps]
        colors = phases
        
        scatter = self.ax2.scatter(positions, [0]*N, s=sizes, c=colors, cmap='hsv', vmin=-180, vmax=180, 
                                   edgecolors='darkblue', linewidth=2, alpha=0.8)
        self.ax2.set_title("Disposición de Antenas en el Arreglo Lineal", fontsize=12, fontweight='bold', pad=10)
        self.ax2.set_xlabel('Posición (λ)', fontsize=10)
        self.ax2.set_ylim(-0.25, 0.25)
        self.ax2.set_yticks([])
        self.ax2.grid(True, axis='x', alpha=0.3)
        
        # Etiquetas informativas
        for i, (pos, amp, phase) in enumerate(zip(positions, amps, phases)):
            self.ax2.text(pos, 0.12, f'A{i+1}', ha='center', va='bottom', fontsize=9, fontweight='bold')
            self.ax2.text(pos, 0.05, f'A:{amp:.2f}', ha='center', va='top', fontsize=8)
            self.ax2.text(pos, -0.05, f'φ:{phase:.0f}°', ha='center', va='bottom', fontsize=8)
        
        # Barra de color
        cbar = self.fig.colorbar(scatter, ax=self.ax2, orientation='horizontal', pad=0.15, shrink=0.8, aspect=30)
        cbar.set_label('Fase (grados)', fontsize=9)
    
    def on_mouse_move(self, event):
        if hasattr(self, 'ax1') and event.inaxes == self.ax1:
            try:
                # Funcionalidad del cursor en 2D
                if self.plot_type.get() == "polar":
                    theta_mouse = event.theta if hasattr(event, 'theta') else np.arctan2(event.ydata, event.xdata)
                else:  # cartesiano
                    if event.xdata is not None and event.ydata is not None:
                        theta_mouse = np.arctan2(event.ydata, event.xdata)
                    else:
                        return
                if 0 <= theta_mouse <= np.pi:
                    U_at_theta = self.total_pattern(np.array([theta_mouse]))
                    D_at_theta = 4 * np.pi * (U_at_theta ** 2) / self.integral_U
                    theta_deg = np.rad2deg(theta_mouse)
                    directivity_val = D_at_theta[0]
                    self.cursor_info.config(text=f"θ: {theta_deg:.1f}° | D: {directivity_val:.2f}")
            except:
                pass

    def update_antenna_widgets(self, *args):
        """Actualizar visibilidad de entradas por antena según N"""
        N = self.N.get()
        for i in range(10):
            if i < N:
                for widget in self.ant_widgets[i]:
                    widget.grid()
            else:
                for widget in self.ant_widgets[i]:
                    widget.grid_remove()

# Demo: Lanzar GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = AntennaArrayGUI(root)
    root.mainloop()