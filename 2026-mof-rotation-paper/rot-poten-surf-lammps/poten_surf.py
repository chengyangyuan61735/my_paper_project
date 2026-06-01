import re
import matplotlib.pyplot as plt
import numpy as np

# 解析LAMMPS输出文件
def parse_lammps_output(filename):
    pattern = r"Rotation angle: (-?\d+) deg, Potential Energy: (\d+\.\d+)"
    angles = []
    energies = []
    
    with open(filename, 'r') as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                angle = float(match.group(1))
                energy = float(match.group(2))
                angles.append(angle)
                energies.append(energy)
    
    return np.array(angles), np.array(energies)

# 绘制势能曲线并保存数据
def plot_and_save_energy_curve(angles, energies, output_file="rotation_energy.png", data_file="energy_data.txt"):
    # 计算相对能量 (以最小能量为0点)
    min_energy = np.min(energies)
    relative_energy = energies - min_energy
    
    # 创建图表
    plt.figure(figsize=(10, 6))
    plt.plot(angles, relative_energy, 'o-', linewidth=2, markersize=8)
    
    # 标记能量最低点
    min_idx = np.argmin(relative_energy)
    plt.annotate(f'Min: {relative_energy[min_idx]:.3f} kcal/mol', 
                 (angles[min_idx], relative_energy[min_idx]),
                 xytext=(10, 20), textcoords='offset points',
                 arrowprops=dict(arrowstyle="->", color='red'))
    
    # 标记能量最高点
    max_idx = np.argmax(relative_energy)
    plt.annotate(f'Max: {relative_energy[max_idx]:.3f} kcal/mol', 
                 (angles[max_idx], relative_energy[max_idx]),
                 xytext=(10, -30), textcoords='offset points',
                 arrowprops=dict(arrowstyle="->", color='blue'))
    
    # 计算势垒高度
    barrier = relative_energy[max_idx] - relative_energy[min_idx]
    plt.title(f'MIL-47 Ligand Rotation Energy Barrier: {barrier:.3f} kcal/mol', fontsize=14)
    
    # 设置坐标轴标签
    plt.xlabel('Rotation Angle (degrees)', fontsize=12)
    plt.ylabel('Relative Energy (kcal/mol)', fontsize=12)
    plt.grid(linestyle='--', alpha=0.7)
    
    # 保存图表
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved energy curve to: {output_file}")
    
    # 保存数据
    with open(data_file, 'w') as f:
        f.write("Angle(deg)  Energy(kcal/mol)  Relative_Energy(kcal/mol)\n")
        for a, e, re in zip(angles, energies, relative_energy):
            f.write(f"{a:>8} {e:>18.6f} {re:>25.6f}\n")
    print(f"Saved energy data to: {data_file}")
    
    plt.show()


if __name__ == "__main__":

    lammps_output_file = "slurm-5721.out"  
    angles, energies = parse_lammps_output(lammps_output_file)
    sorted_indices = np.argsort(angles)
    sorted_angles = angles[sorted_indices]
    sorted_energies = energies[sorted_indices]
    plot_and_save_energy_curve(sorted_angles, sorted_energies)