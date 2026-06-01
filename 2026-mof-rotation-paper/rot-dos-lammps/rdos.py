import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from pathlib import Path
import re

def read_omega_file(filename):
    """读取 omega.out 文件，返回角速度数据"""
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    timestep_data = []
    while i < len(lines):
        line = lines[i].strip()
        
        
        if line.startswith('#') or not line:
            i += 1
            continue
        
        # 处理时间步行 (如 "500000 97")
        if re.match(r'^\d+\s+\d+$', line):

            if timestep_data:
                data.append(np.array(timestep_data))
                timestep_data = []
            
            # 读取时间步和行数
            timestep, n_rows = map(int, line.split())
            
            # 读取接下来的n_rows行
            for j in range(n_rows):
                i += 1
                if i >= len(lines):
                    break
                row_line = lines[i].strip()
                if row_line:
                    parts = row_line.split()
                    if len(parts) >= 4:
                        # 只保留角速度数据，忽略行号
                        omega_vals = list(map(float, parts[1:4]))
                        timestep_data.append(omega_vals)
        i += 1
    
    # 添加最后一个时间步的数据
    if timestep_data:
        data.append(np.array(timestep_data))
    
    return np.array(data)

def compute_acf(omega_data, direction, Nc=1000):
  
    # 提取指定方向的数据 (Nf, n_chunks)
    dir_index = ['x', 'y', 'z'].index(direction.lower())
    omega_dir = omega_data[:, :, dir_index]
    
    Nf, n_chunks = omega_dir.shape
    M = Nf - Nc  # 时间原点数量
    
    
    acf = np.zeros((n_chunks, Nc))
    for nc in range(Nc):
        # 使用矩阵运算替代三重循环
        acf[:, nc] = np.mean(omega_dir[:M] * omega_dir[nc:nc+M], axis=0)
    
    # 对所有chunk取平均
    acf_avg = np.mean(acf, axis=0)
    return acf_avg

def compute_dos(acf, dt=0.001, Nc=1000, omega_max=1000):
    """计算旋转态密度(PDOS)"""
    # 归一化自相关函数
    vacf = acf / acf[0]
    
    # 加窗函数
    window = (np.cos(np.pi * np.arange(Nc) / Nc) + 1) * 0.5
    vacf_windowed = vacf * window
    
    
    vacf_sym = vacf_windowed.copy()
    vacf_sym[1:] *= 2
    
    
    omega_vals = np.arange(0, omega_max + 1)
    
    # FFT
    t = np.arange(Nc) * dt
    cos_terms = np.cos(np.outer(omega_vals, t))
    pdos = dt * np.dot(cos_terms, vacf_sym)
    
    # 转换为频率 (THz)
    nu_vals = omega_vals / (2 * np.pi)
    return nu_vals, pdos

def process_single_file(filename, direction='z', output_dir='results'):
    """处理单个文件并保存结果"""
    
    os.makedirs(output_dir, exist_ok=True)
    base_name = Path(filename).stem
    
    print(f"处理文件: {filename} ({direction}方向)")
    
    
    omega_data = read_omega_file(filename)
    
    
    acf = compute_acf(omega_data, direction)
    
    
    nu_vals, pdos = compute_dos(acf)
    
    # 保存结果
    acf_filename = f"{output_dir}/ACF_{base_name}_{direction}.txt"
    np.savetxt(acf_filename, np.column_stack((np.arange(len(acf)) * 0.001, acf)), 
               header='correlation_time(ps) ACF', fmt='%.6f')
    
    dos_filename = f"{output_dir}/DOS_{base_name}_{direction}.txt"
    np.savetxt(dos_filename, np.column_stack((nu_vals, pdos)), 
               header='Frequency(THz) PDOS', fmt='%.6f')
    
    # 绘图
    plt.figure(figsize=(12, 10))
    
    plt.subplot(2, 1, 1)
    plt.plot(np.arange(len(acf)) * 0.001, acf, 'b-', linewidth=2)
    plt.xlabel('Correlation Time (ps)', fontsize=12)
    plt.ylabel(f'ACF of Angular Velocity ({direction.upper()})', fontsize=12)
    plt.title(f'Autocorrelation Function - {base_name}', fontsize=14)
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.plot(nu_vals, pdos, 'r-', linewidth=2)
    plt.xlabel('Frequency (THz)', fontsize=12)
    plt.ylabel('Rotational DOS', fontsize=12)
    plt.title(f'Rotational Density of States - {base_name}', fontsize=14)
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/Results_{base_name}_{direction}.png")
    plt.close()
    
    return nu_vals, pdos

def process_files(file_pattern, direction='z', output_dir='results'):
    """处理匹配文件模式的所有文件并计算平均PDOS"""
    
    file_list = sorted(glob.glob(file_pattern))
    
    if not file_list:
        print(f"未找到匹配的文件: {file_pattern}")
        return
    
    print(f"找到 {len(file_list)} 个匹配文件:")
    for f in file_list:
        print(f"  - {f}")
    
    all_pdos = []
    base_names = []
    nu_vals = None
    
    
    for filename in file_list:
        base_name = Path(filename).stem
        base_names.append(base_name)
        
        # 处理单个文件
        try:
            current_nu_vals, pdos = process_single_file(filename, direction, output_dir)
            
            # 确保所有文件的频率范围一致
            if nu_vals is None:
                nu_vals = current_nu_vals
            elif not np.array_equal(nu_vals, current_nu_vals):
                print(f"警告: 文件 {filename} 的频率范围与其他文件不一致")
            
            all_pdos.append(pdos)
        except Exception as e:
            print(f"处理文件 {filename} 时出错: {str(e)}")
            continue
    
    if not all_pdos:
        print("没有成功处理任何文件")
        return
    
    # 计算平均PDOS
    avg_pdos = np.mean(all_pdos, axis=0)
    
    # 保存平均结果
    avg_filename = f"{output_dir}/AVG_DOS_{direction}.txt"
    np.savetxt(avg_filename, np.column_stack((nu_vals, avg_pdos)), 
               header='Frequency(THz) PDOS', fmt='%.6f')
    
    # 绘制所有曲线和平均曲线
    plt.figure(figsize=(10, 6))
    for i, pdos in enumerate(all_pdos):
        plt.plot(nu_vals, pdos, alpha=0.5, label=base_names[i])
    
    plt.plot(nu_vals, avg_pdos, 'k-', linewidth=3, label='Average')
    plt.xlabel('Frequency (THz)', fontsize=12)
    plt.ylabel('Rotational DOS', fontsize=12)
    plt.title(f'Rotational Density of States ({direction.upper()} Direction)', fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/All_DOS_{direction}.png")
    plt.close()
    
    return nu_vals, avg_pdos


if __name__ == "__main__":
    # 处理单个文件
    # process_single_file('omega1.out', direction='z')
    
    # 处理多个文件 (使用通配符)
    process_files('omega*.out', direction='x')
