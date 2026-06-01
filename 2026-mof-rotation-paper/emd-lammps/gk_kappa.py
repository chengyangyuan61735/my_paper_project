import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as kB
import glob
import os

class ThermalConductivityAnalyzer:
    def __init__(self, config):
        self.config = config
        self.results = {}
        
    def _read_log_file(self, file_path):
        """read temp and volume data from log file"""
        try:
            df = pd.read_csv(file_path, 
                             skiprows=self.config['head'],
                             nrows=self.config['ll'],
                             delim_whitespace=True,
                             header=None,
                             names=['Step', 'Temp', 'Press', 'Volume', 'KinEng', 
                                    'PotEng', 'TotEng', 'Lx', 'Ly', 'Lz', 'Density', 'CPU'])
            
            T_avg = df['Temp'].mean()
            V0 = df['Volume'].iloc[0]
            return T_avg, V0
        except Exception as e:
            print(f"Error reading {file_path}: {str(e)}")
            return None, None

    def _read_hafc_file(self, file_path):
        try:
            df = pd.read_csv(file_path, 
                             skiprows=3, 
                             delim_whitespace=True,
                             header=None,
                             names=['index', 'timedelta', 'ncount', 'x', 'y', 'z'])
            return df[['x', 'y', 'z']].tail(self.config['n']).values.T
        except Exception as e:
            print(f"Error reading {file_path}: {str(e)}")
            return None

    def _integrate_heat_flux(self, j_components):
        n = j_components.shape[1]
        k_components = np.zeros_like(j_components)
        
        for i in range(3):
            j = j_components[i]
            cumsum = np.cumsum(j)
            k_components[i] = cumsum - 0.5*j[0] - 0.5*j
            
        return k_components

    def calculate_scale(self, log_file):
        T_avg, V0 = self._read_log_file(log_file)
        if T_avg is None or V0 is None:
            return None
            
        dt = self.config['dt']
        convert = self.config['convert']
        scale = dt * convert / (kB * T_avg**2 * V0)
        return scale

    def process_case(self, hafc_file, scale, direction='average'):
        j_xyz = self._read_hafc_file(hafc_file)
        if j_xyz is None:
            return None
            
        k_xyz = self._integrate_heat_flux(j_xyz)
        
        if direction == 'x':
            kappa = scale * k_xyz[0]
        elif direction == 'y':
            kappa = scale * k_xyz[1]
        elif direction == 'z':
            kappa = scale * k_xyz[2]
        else:
            kappa = scale * np.mean(k_xyz, axis=0)
        
        return kappa

    def visualize_results(self, kappa_curves):
        plt.figure(figsize=(12, 7))
        t = np.arange(1, self.config['n']+1)
        num_files = len(kappa_curves)
        
        if num_files <= 8:
            colors = ['red', 'blue', 'green', 'cyan', 'magenta', 'yellow', 'orange', 'purple']
        else:
            colors = plt.cm.rainbow(np.linspace(0, 1, num_files))
        
        for i, kappa in enumerate(kappa_curves):
            if kappa is not None:
                label = f'Seed {i+1}' if num_files <= 20 else None
                plt.plot(t, kappa, color=colors[i % len(colors)], alpha=0.7, 
                         label=label)
        
        # plot average 
        valid_curves = [k for k in kappa_curves if k is not None]
        if valid_curves:
            kappa_avg = np.mean(valid_curves, axis=0)
            plt.plot(t, kappa_avg, 'k-', linewidth=2.5, label='Average')
        
        plt.title(f'Thermal Conductivity ({num_files} seeds)')
        plt.xlabel('Time (fs)')
        plt.ylabel('Thermal Conductivity (W/mK)')
        
        if num_files <= 20:
            plt.legend()
        else:
            plt.legend(['Average'])
            
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        # convergence interval selected by the user
        print("Please select two points on the plot to define the convergence region")
        points = plt.ginput(2, timeout=-1)
        plt.close()
        
        return points

    def calculate_converged_kappa(self, kappa_curves, points):
        start_idx = max(0, int(round(points[0][0])))
        end_idx = min(self.config['n'], int(round(points[1][0])))
        
        kappa_values = []
        for i, kappa in enumerate(kappa_curves):
            if kappa is not None:
                kappa_conv = np.mean(kappa[start_idx:end_idx])
                kappa_values.append(kappa_conv)
        
        kappa_avg = np.mean(kappa_values)
        kappa_std = np.std(kappa_values, ddof=1) / (len(kappa_values)**0.5)
        
        return kappa_values, kappa_avg, kappa_std

    def save_results(self, kappa_values, kappa_avg, kappa_std):
        num_files = len(kappa_values)
        
        with open('kappa_z.txt', 'w') as f:
            header = " ".join([f"kappa{i+1}" for i in range(num_files)]) 
            header += " kappa_ave kappa_error_bar"
            f.write(header + "\n")
            
            values_str = " ".join(f"{k:.6g}" if k is not None else "NaN" for k in kappa_values)
            f.write(f"{values_str} {kappa_avg:.6g} {kappa_std:.6g}")

    def find_matching_files(self, pattern):
        files = sorted(glob.glob(pattern))
        if not files:
            print(f"No files found matching pattern: {pattern}")
            return []
            
        print(f"Found {len(files)} files: {files[:3]}...")  # 显示前3个文件
        return files

# preset parameters
config = {
    'n': 20001,       # heat flux correlation steps
    'dt': 5,          # timestep, dt * output_step
    'convert': 4.8270e-16,  # unit convert
    'll': 1000,       # log read line
    'head': 1203,     # omit line
}

analyzer = ThermalConductivityAnalyzer(config)

# user define file number and direction
use_default = input("Use default file patterns? (y/n): ").strip().lower()
if use_default == 'y':
    log_pattern = "log.*"
    hafc_pattern = "hafc*.txt"
    num_files = int(input("Number of files to process: ").strip() or 8)
else:
    log_pattern = input("Enter log file pattern (e.g. 'log.*' or 'log.[1-50]'): ").strip()
    hafc_pattern = input("Enter HAFC file pattern (e.g. 'hafc*.txt'): ").strip()
    num_files = None  

log_files = analyzer.find_matching_files(log_pattern)
hafc_files = analyzer.find_matching_files(hafc_pattern)

if num_files is not None:
    log_files = log_files[:num_files]
    hafc_files = hafc_files[:num_files]

direction = input("Select direction (x, y, z, average): ").strip().lower()

scales = [analyzer.calculate_scale(log_file) for log_file in log_files]

kappa_curves = []
for hafc_file, scale in zip(hafc_files, scales):
    if scale is not None:
        kappa = analyzer.process_case(hafc_file, scale, direction)
        kappa_curves.append(kappa)
    else:
        kappa_curves.append(None)

points = analyzer.visualize_results(kappa_curves)

kappa_values, kappa_avg, kappa_std = analyzer.calculate_converged_kappa(kappa_curves, points)

# save results
analyzer.save_results(kappa_values, kappa_avg, kappa_std)

# print
print(f"\nProcessed {len(log_files)} files")
print(f"Selected convergence region: {points[0][0]:.0f} - {points[1][0]:.0f} fs")
print(f"Direction: {direction.upper() if direction != 'average' else 'Average (x+y+z)/3'}")
print(f"Average thermal conductivity: {kappa_avg:.4f} ± {kappa_std:.4f} W/mK")
print(f"Results saved to kappa.txt")
