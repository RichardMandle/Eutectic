import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import importlib

'''
Simple python script for calculating phase diagrams and eutectic points.

Functions:
save_data                                - saves data to .csv file (warning, files can be huge)
build_ternary_axis                       - called during plotting; makes the ternary axis system used

print_composition                       - prints the composition of the eutectic point (default) or ...
make_phase_diagram                      - Phase diagram engine - calculates the melting point at each concentration
plot_phase_diagram                      - Plotting engine - makes plots for 2- or 3- component systems.

'''


def save_data(concentrations, melt_points, clear_points, filename='EutecticData.csv'):
    '''
    Save data to CSV file
    headers tell you whats what.
    '''
    base, ext = os.path.splitext(filename)
    suffix = 1
    while os.path.exists(filename):
        filename = f"{base} ({suffix}){ext}"
        suffix += 1

    headers = [f"conc {i+1}" for i in range(concentrations.shape[1])]
    headers.extend(["melt_points", "clear_points"])

    data = np.column_stack((concentrations, melt_points, clear_points))
    np.savetxt(filename, data, delimiter=",", header=",".join(headers), fmt="%.5f")
    print(f'Saved data to: {filename}')
    return

def print_composition(concentrations, melt_points, clear_points, quiet=True):
    
    # Prints eutectic composition and phase transition data.
    
    idx = np.argmin(melt_points)
    melting_point = np.round(melt_points[idx], 2)
    clearing_point = np.round(clear_points[idx], 2)
    enantiotropic_range = np.round(clearing_point - melting_point, 2)

    if not quiet:
        print('\nEutectic composition:')
        for i, conc in enumerate(concentrations[idx]):
            comp_percentage = np.round(100 * conc, 3)
            print(f"{comp_percentage}% Component {chr(ord('A') + i)}")
        print(f'Melting point: {melting_point} °C')
        print(f'Clearing point: {clearing_point} °C')
        print(f'Enantiotropic range: {enantiotropic_range} °C')
    return melting_point, clearing_point, enantiotropic_range

def make_phase_diagram(melting_points, enthalpy_fusion, clearing_points,
                       search_size=25000, plotting=True, quat_plot_style=0):
    
    #Samples mixture compositions and calculates eutectic and clearing points.
    
    melting_points = np.array(melting_points)
    enthalpy_fusion = np.array(enthalpy_fusion)
    clearing_points = np.array(clearing_points)
    R = 8.314  # J/(mol·K)

    n_components = len(melting_points)
    random_compositions = np.random.rand(search_size, n_components)
    concentrations = random_compositions / random_compositions.sum(axis=1, keepdims=True)

    if np.any(concentrations <= 0):
        raise ValueError("Zero or negative concentrations encountered.")

    T_melt_K = np.tile(melting_points + 273.15, (search_size, 1))
    dH_fusion_J = np.tile(enthalpy_fusion * 1000, (search_size, 1))
    T_clear_C = np.tile(clearing_points, (search_size, 1))

    melt_points = np.max(dH_fusion_J / ((dH_fusion_J / T_melt_K) - R * np.log(concentrations)), axis=1) - 273.15
    clear_points = np.sum(T_clear_C * concentrations, axis=1)
    
    print_composition(concentrations, melt_points, clear_points, quiet = False)
    
    if plotting:
        plot_phase_diagram(concentrations, melt_points, clear_points, quat_plot_style)

    return concentrations, melt_points, clear_points

def plot_phase_diagram(concentrations, melt_points, clear_points, quat_plot_style = 0):
    
    # plotting for 2, 3, or 4 component systems...
    n_components = concentrations.shape[1]

    if n_components == 2:
        plt.figure()
        plt.scatter(concentrations[:, 1], melt_points, label='Melting Point (°C)', alpha=0.6)
        plt.scatter(concentrations[:, 1], clear_points, label='Clearing Point (°C)', alpha=0.6)
        plt.xlabel('Mole Fraction of Component B')
        plt.ylabel('Temperature (°C)')
        plt.title('Binary Phase Diagram')
        plt.legend()
        plt.grid(True)
        plt.xlim([0, 1])

        idx = np.argmin(melt_points)
        plt.axvline(concentrations[idx, 1], color='k', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.show()

    elif n_components == 3:
        tri_plot = importlib.import_module("tri_plot")
        tri_plot.mpl_tri_plot(concentrations, melt_points, clear_points)

    elif n_components == 4:
        if quat_plot_style == 0:
            print("Plotting skipped for quaternary system.")
        else:
            print("Loading quaternary plotting module, please wait...")
            quat_plot = importlib.import_module("quat_plot")
            print(f"Attempting pyramidal plot for {n_components} components")

            style_map = {
                1: quat_plot.mpl_quat_plot,
                2: quat_plot.mvi_quat_scatter,
                3: quat_plot.mvi_quat_isosurf
            }

            if quat_plot_style in style_map:
                style_map[quat_plot_style](concentrations, melt_points, clear_points)
            else:
                print("Unknown plotting style for quaternary plot.")

    else:
        print(f"No plotting for {n_components}-component system.")


# function to get user input
def get_input(prompt):
    return input(prompt).strip()

# function to parse input strings into lists of floats
def parse_input(input_str):
    return [float(x.strip()) for x in input_str.split(',') if x.strip()]

def main():
    '''
    if called in this way, prompt the user to give some inputs we can use
    '''
    print("\n--- Eutectic Diagram Generator ---\n")

    melting_points = parse_input(get_input("Enter melting points (°C), comma-separated: "))
    enthalpies = parse_input(get_input("Enter enthalpies of fusion (kJ/mol), comma-separated: "))

    if len(melting_points) != len(enthalpies):
        print("ERROR: Mismatched melting points and enthalpies!")
        return

    clearing_points = parse_input(get_input("Enter clearing points (°C), comma-separated: "))
    if len(melting_points) != len(clearing_points):
        print("ERROR: Mismatched melting points and clearing points!")
        return

    plotting_style = 0
    if len(melting_points) == 4:
        plotting_style = int(get_input("Select quaternary plot style:\n 0 = None\n 1 = Matplotlib\n 2 = Mayavi Scatter\n 3 = Mayavi Isosurface\n> "))

    concentrations, melt_points, clear_points = make_phase_diagram(
        melting_points, enthalpies, clearing_points, quat_plot_style=plotting_style
    )
   

    if get_input("\nSave results to CSV? (y/n): ").lower() == 'y':
        save_data(concentrations, melt_points, clear_points)

if __name__ == "__main__":
    main()
