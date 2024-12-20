# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import os, re

# Function to extract relevant data from an AMS output file
def ext_xtb(name='test.out'):
    sp = 0  # Flag to indicate single-point calculation
    g_e, freq, Nmods, geo, geo_sp, coord, counter, bonds = [[], []], [], [], [], [], "", 0, []

    # Read the file contents
    with open(name) as f:
        data = f.readlines()

    # Split the lines into words for easy processing
    reci = [[j for j in i[:-1].split()] for i in data]

    # Find the number of atoms (NA) from the output
    for i in range(len(data)):
        if 'Total System Charge' in data[i]:
            NA = int(data[i-2].split()[0])  # Number of atoms appears 2 lines above
            break

    # Extract the last geometry (in Bohr units) and corresponding single-point geometry (in Angstroms)
    for i in range(len(data)-1, 0, -1):
        if "Index Symbol   x (angstrom)   y (angstrom)   z (angstrom)" in data[i]:
            for k in data[i+1:i+NA+1]:
                # Convert coordinates to Bohr (multiplying by 1.88973)
                geo += [" ".join([str(float(k[12:].split()[j]) * 1.88973) if j > 0 else k[12:].split()[j] for j in range(4)]) + '\n']
                geo_sp += [" ".join([k[12:].split()[j] for j in range(4)]) + '\n']
            break

    # Extract optimization steps or transition state search data
    for i in range(len(data)): 
        if "*  GEOMETRY OPTIMIZATION" in data[i] or '*  TRANSITION STATE SEARCH  *' in data[i]:                            
            for j in range(len(data)):
                if "Geometry Convergence after Step" in data[j]:
                    counter += 1
                    # Extract energy (in kcal/mol) and gradient values
                    g_e[1] += [float(reci[j+2][2]) * 627.51]  # Energy in kcal/mol
                    g_e[0] += [float(reci[j+4][3])]  # Gradient
                    # Create coordinate block for visualization
                    coord += f'{NA}\n{counter}_{g_e[0][-1]}\n'
                    for k in data[j+10:j+24+NA]:
                        if re.search(r'-?\d\.\d{6}\d*\s+-?\d\.\d{6}\d*', k):
                            coord += k[12:]
            g_e = np.array(g_e)
            g_e[1, :] = g_e[1, :] - np.min(g_e[1, :])  # Normalize energies
            
            # Additional processing for TS searches
            if '*  TRANSITION STATE SEARCH  *' in data[i]:
                data_ = [i[:-1] for i in data if 'Distance' in i][0].split()[2:4]
                if int(data_[0]) > int(data_[1]):  # Ensure correct atom order
                    data_ = [data_[1], data_[0]]
                with open('test.out') as f:
                    bonds = np.array([float(i[:-1].split()[6]) for i in f if f"bnd{data_[0]:>8}{data_[1]:>8}" in i])
        elif '* SINGLE POINT CALCULATION *' in data[i]:
            sp = 1  # Single-point calculation flag
            with open('last.angs', 'w') as f:
                for red in geo_sp:
                    print(red[:-1], file=f)
            with open('last.xyz', 'w') as f:
                print(len(geo_sp), '\n', file=f)
                for red in geo_sp:
                    print(red[:-1], file=f)

        # Extract frequency information and save in Molden format
        if "Normal Mode Frequencies" in data[i]:
            freq = [j.split()[1] + '\n' for j in data[i+4:i+3*NA-2]]
            for k in range(len(freq)):
                Nmods += [f'vibration    {k+1}\n']
                Nmods += [l[16:] for l in data[i+3*NA+8+k*(NA+3):i+4*NA+11+k*(NA+3)] if re.search(r'-?\d\.\d{6}', l)]
            molden = ['[Molden Format]\n', '[Atoms] AU\n', '[FREQ]\n', '[FR-COORD]\n', '[FR-NORM-COORD]\n']
            molden_format = molden[:2] + geo + [molden[2]] + freq + [molden[3]] + geo + [molden[4]] + Nmods
            with open('frq.molden', 'w') as f:
                for line in molden_format:
                    f.write(line)
            break

        # Handle cases where geometry did not converge
        if "Last geometry (not converged)" in data[i]:
            print(f'Geometry not converged! minG {round(g_e[0,:][np.argmin(g_e[0,:])], 4)}, at step {np.argmin(g_e[0,:])}, while last step is {round(g_e[0,:][-1], 4)}\n{os.getcwd()}')
    
    # Save optimization trajectory
    with open('opt.xyz', 'w') as f:
        f.write(coord)
    
    # Save the last geometry
    if not sp:
        koordinate = coord.split('\n')
        with open('last.angs', 'w') as f:
            for line in koordinate[-(NA+1):-1]:
                print(line, file=f)
        with open('last.xyz', 'w') as f:
            for line in koordinate[-(NA+3):-1]:
                print(line, file=f)

    return g_e, coord, bonds, freq, geo, Nmods

# Main script execution
with open('test.out') as f:
    output = f.read()
g_e, C, bonds, freq, geo, Nmods = ext_xtb()

# Plotting optimization or TS search results
if "*  GEOMETRY OPTIMIZATION" in output or '*  TRANSITION STATE SEARCH  *' in output or '* Potential Energy Surface Scan *' in output:
    fig, ax1 = plt.subplots() 
    ax1.set_xlabel('step') 
    ax1.set_ylabel('E, kcal/mol', color='red', rotation="horizontal") 
    ax1.plot(np.arange(1, g_e[1, :].size+1), g_e[1, :], '|-', color=(1, 0, 0, 0.5))
    plt.grid()
    
    # Add gradient to the plot
    if "*  GEOMETRY OPTIMIZATION" in output or '*  TRANSITION STATE SEARCH  *' in output:
        ax2 = ax1.twinx() 
        ax2.set_ylabel('grad_max', color=(0, 0, 1, 0.6), loc="top", rotation="horizontal") 
        ax2.plot(np.arange(1, g_e[1, :].size+1), g_e[0, :], '--', color=(0, 0, 1, 0.3)) 
        ax2.tick_params(axis='y', colors=(0, 0, 1, 0.4), direction='in', pad=-22)
        
        # Add reaction coordinate for TS search
        if '*  TRANSITION STATE SEARCH  *' in output:
            ax3 = ax1.twinx() 
            ax3.plot(np.arange(1, g_e[1, :].size+1), bonds[1:], '-', color=(1, 0.5, 0.15, 0.3)) 
            ax3.set_ylabel('TS coordinate', color=(1, 0.5, 0.15, 0.7), rotation="horizontal") 
            ax3.tick_params(axis='y', colors=(1, 0.5, 0.15, 0.5))
    
    plt.savefig('opt.png')  # Save the plot
