import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import ticker
from matplotlib.patches import Circle, PathPatch
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

import string
import os
import imageio
import sys

'''
Simple python script for calculating phase diagrams and eutectic points.

Functions:
save_data                                - saves data to .csv file (warning, files can be huge)
build_ternary_axis                       - called during plotting; makes the ternary axis system used

plot_pyramidal_ax                       - generates a triangular pyramidal axis for 4-component plotting
label_points                            - labels the pyramidal axis
get_cartesian_array_from_barycentric    - gets a cartesian array from the barycentre of the pyramidal plot
plot_3d_tern                            - plots a 3D scatter plot within a pyramidal axis
get_random_subset                       - generates a random subset of compositions which includes the eutectic

print_composition                       - prints the composition of the eutectic point (default) or ...
make_phase_diagram                      - Phase diagram engine - calculates the melting point at each concentration
plot_phase_diagram                      - Plotting engine - makes plots for 2- or 3- component systems.

'''

def save_data(concentrations, melt_points, clear_points, filename='EutecticData.csv'):
    '''
    Save data to CSV file
    headers tell you whats what.
    '''
    if os.path.exists(filename): # check if filename exists
        # If it does, add a suffix to the filename to make it unique
        suffix = 1
        while True:
            new_filename = "{} ({}){}".format(*os.path.splitext(filename) + (suffix,))
            if not os.path.exists(new_filename):
                filename = new_filename
                break
            suffix += 1

    headers = ["conc {}".format(i) for i in range(1,concentrations.shape[1]+1)]
    headers.extend(["melt_points", "clear_points"])

    data = np.column_stack((concentrations,melt_points,clear_points))     # Stack the arrays horizontally

    np.savetxt(filename, data, delimiter=",", header=",".join(headers), fmt="%s")     # Save the data to a CSV file

    print('saved data to: ' + filename)

    return


def build_ternary_axis(divisions,z_values):
    # Builds the triangular axis system used for ternary plots

    d1 = np.cos(np.pi/3)
    d2 = np.sin(np.pi/3)

    l = np.linspace(0,1,divisions)
    zmax = max(z_values)
    for i in range (0,np.size(l)):
            plt.plot([l[i]*d1, 1-l[i]*d1], [l[i]*d2, l[i]*d2],color='black',linewidth='0.5') # horizontal lines
            plt.plot([l[i], l[i]+(1-l[i])*d1],[0, (1-l[i])*d2], color='black',linewidth='0.5') # upper right->lower left lines)
            plt.plot([(1-l[i])*d1, 1-l[i]],[(1-l[i])*d2, 0], color='black',linewidth='0.5') # upper left->lower right lines))

            plt.text(l[i],-0.1, np.array2string(np.around(l[i],decimals=3))) #x axis labels
            plt.text((1-l[i])*np.cos(np.pi/3)+0.525, l[i]*np.sin(np.pi/3), np.array2string(np.around(l[i],decimals=3))) #right axis labels
            plt.text(0.5-l[i]*np.cos(np.pi/3)-0.15, np.sin(np.pi/3)*(1-l[i]), np.array2string(np.around(l[i],decimals=3))) #left axis labels

    plt.text(0.20,-0.155,'Component A / mol %',rotation=0)
    plt.text(-0.125,0.65,'Component B / mol %',rotation=60)
    plt.text(0.775,0.65,'Component C / mol %',rotation=-60)


def plot_pyramidal_ax(axis):
    """
    Plot the tetrahedral outline in 3D.

    The function plots the tetrahedral outline defined by the vertices
    using lines connecting the vertices.

    No inputs, no returns
    """

    verts = [[0, 0, 0],
             [1, 0, 0],
             [0.5, np.sqrt(3) / 2, 0],
             [0.5, 0.28867513, 0.81649658]]
    lines = combinations(verts, 2)
    for x in lines:
        line = np.transpose(np.array(x))
        axis.plot3D(line[0], line[1], line[2], c='0')


def label_points(axis):
    """
    Label the vertices of the tetrahedron in the plot.

    The function labels each vertex of the tetrahedron using the Barycentric
    coordinates and assigns a corresponding label ('a', 'b', 'c', 'd').

    No inputs, no returns
    """

    a = np.array([1, 0, 0, 0])
    b = np.array([0, 1, 0, 0])
    c = np.array([0, 0, 1, 0])
    d = np.array([0, 0, 0, 1])
    labels = ['a', 'b', 'c', 'd']
    cartesian_points = get_cartesian_array_from_barycentric(np.array([a, b, c, d]))
    for point, label in zip(cartesian_points, labels):
        if 'a' in label:
            axis.text(point[0], point[1] - 0.075, point[2], label, size=16)
        elif 'b' in label:
            axis.text(point[0] + 0.02, point[1] - 0.02, point[2], label, size=16)
        else:
            axis.text(point[0], point[1], point[2], label, size=16)


def get_cartesian_array_from_barycentric(b):
    """
    Convert Barycentric coordinates to Cartesian coordinates.

    The function takes an array of Barycentric coordinates and converts
    them to Cartesian coordinates based on the predefined vertices.

    Parameters:
    - b (ndarray): Array of Barycentric coordinates.

    Returns:
    - t_array (ndarray): Array of Cartesian coordinates.
    """
    verts = [[0, 0, 0],
             [1, 0, 0],
             [0.5, np.sqrt(3) / 2, 0],
             [0.5, 0.28867513, 0.81649658]]
    t = np.transpose(np.array(verts))
    t_array = np.array([t.dot(x) for x in b])
    return t_array


def plot_3d_tern(axis,concs, melting_points, alpha_array):
    """
    Plot 3D scatter points in the ternary coordinate system.

    The function plots the scatter points based on the given concentration
    values (concs) in the ternary coordinate system. The color of the points
    is determined by the melting points (melting_points).

    Parameters:
    - concs (ndarray): Array of concentration values.
    - melting_points (ndarray): Array of melting points.
    - alpha_array (ndarray): list of alpha values for points

    Returns:
    None
    """

    cartesian_points = get_cartesian_array_from_barycentric(concs)
    scatter = axis.scatter(cartesian_points[:-1, 0], cartesian_points[:-1, 1], cartesian_points[:-1, 2],
               c=melting_points[:-1], cmap='magma', s=100, alpha=alpha_array)

    #do the eutectic point slightly differently.
    axis.scatter(cartesian_points[-1, 0], cartesian_points[-1, 1], cartesian_points[-1, 2],
               c=melting_points[-1], cmap='magma', s=100, alpha=1)

    return scatter


def get_random_subset(concentrations, melt_points, clear_points, n):
    """
    Obtain a random subset of points and find the eutectic values.

    The function selects a random subset of points from the given
    concentration array (concentrations) and corresponding melting points array (melt_points),
    based on the specified subsampling interval (n). It also identifies the
    lowest value in the subset of melting points and its corresponding
    concentration values.

    random choice will always contain the eutectic

    Parameters:
    - concentrations (ndarray): Array of concentration values.
    - melt_points (ndarray): Array of melting points.
    - n (int): Subsampling interval.

    Returns:
    - random_concs (ndarray): Random subset of concentration values.
    - random_melts (ndarray): Random subset of melting points.

    """

    random_indices = np.append(np.random.choice(range(concentrations.shape[0]), size=concentrations.shape[0] // n, replace=False),
                               np.argmin(melt_points))
    random_concs = np.squeeze(concentrations[random_indices,:])
    random_melts = np.squeeze(melt_points[random_indices])
    random_clear = np.squeeze(clear_points[random_indices])

    return random_concs, random_melts, random_clear

def print_composition(concentrations, melt_points, clear_points, quiet = False):
    '''
    Simply prints the eutectic composition, i.e. the lowest melt, from an array of concentrations and corresponding melting/clearing points

    Inputs:
    concentrations - an array of concentrations for each mixture
    melt_points - a list of melting points for each mixture (in °C)
    clear_points - a list of clearing points for each mixture (in °C)

    returns:
    some text
    '''
    idx = np.argmin(melt_points)
    melting_point = np.round(melt_points[idx], decimals=2)
    clearing_point = np.round(clear_points[idx], decimals=2)
    enantiotropic_range = np.round(clearing_point - melting_point, decimals=2)

    if not quiet:
        print ('Eutectic composition:')
        for n in range (0,np.shape(concentrations)[1]):
            comp_percentage = np.round(100 * concentrations[idx, n], decimals=3)
            print (f"{comp_percentage}% Component {chr(ord('A')+n)}")

        print (f'Melting point: {melting_point} °C')
        print (f'Clearing point: {clearing_point} °C')
        print (f'Enantiotropic range: {enantiotropic_range} °C')

    return melting_point, clearing_point, enantiotropic_range

def make_phase_diagram(melt_points, enthalpy_fusion, clear_points, search_size=250000, plotting=True):
    '''
    Takes a list of melting points, enthalpies of fusion, molecular weights
    searches a number of random concentrations (set by search_size) and calculates
    the melting point and clearing point of each.

    Works well for non-LCs & nematics, less well for smectics, untested for everything else.

    Basic idea is taken from https://doi.org/10.1039/C39740000098 but implemented differently;
    rather than seek the absolute eutectic we just sample a huge number of different concentrations.

    Inputs:
    melt_points      - the melting point of each component in °C
    enthalpy_fusion       - the enthalpy of fusion of each component (in kJ mol-1)
    clear_points          - the clearing point of each component in °C
    search_size (250000)  - the number of concentrations to search
    Plotting (True)       - whether or not to make a plot; this does slow it down a bit.

    Returns:
    concentrations - the concentrations used; a NxM array, where N is the number of mixtures and M is the number of components
    melt_points  - the melting point of each different concentration
    clear_points - the clearing point of each different concentration
    '''

    temperatures = np.random.rand(search_size, np.size(melt_points))
    concentrations = temperatures / np.transpose(np.tile(np.sum(temperatures, 1), (np.size(melt_points), 1)))

    M = np.double(np.tile(melt_points + 273.15, (search_size, 1)))            # get tile of melting points in K
    E = np.double(np.tile(enthalpy_fusion * 1000, (search_size, 1)))          # ... enthalpies in J/mol
    C = np.double(np.tile(clear_points, (search_size, 1)))                    # ... clearing points in °C

    R = 8.314  # Gas constant in J/(mol·K)

    # Calculate the melting points of the mixtures
    melt_points = np.max(E / ((E / M) - R * np.log(concentrations)), axis=1) - 273.15
    # Calculate the clearing points of the mixtures
    clear_points = np.sum(C * concentrations, axis=1)

    if plotting:
        plot_phase_diagram(concentrations, melt_points, clear_points)

    return concentrations, melt_points, clear_points


def plot_phase_diagram(concentrations, melt_points, clear_points):
    '''
    Plotting engine. Will attempt to make a plot; if we have 2 components it'll be a line graph, if 3 a ternary one.

    Inputs:
    concentrations - array of concentrations of each component in a given mixture
    melt_points  - array of melting points of each mixture
    clear_points - array of clearing points of each mixture

    Returns:
    A plot for 2- or 3- component mixtures.
    '''

    if np.shape(concentrations)[1] == 2:
        fig = plt.figure()
        plt.scatter(concentrations[:,1], melt_points, label='Melting Point / °C')
        plt.scatter(concentrations[:,1], clear_points, label='Clearing Point / °C')
        plt.xlabel('Concentration of B / mol %')
        plt.ylabel('Temperature / °C')
        plt.legend(loc='upper right')
        plt.xlim([0,1])

        # Find the index of the eutectic point
        idx = np.argmin(melt_points)
        # Extract the eutectic composition and corresponding temperatures
        x_eutectic = concentrations[idx, 1]
        y_melting = melt_points[idx]
        y_clearing = clear_points[idx]

        # Prepare data for plotting
        X = [x_eutectic, x_eutectic]
        Y = [y_melting, y_clearing]
        plt.plot(X, Y, color='black', linewidth='0.5', linestyle='dashed')
        plt.show()


    elif np.shape(concentrations)[1] == 3: # if we've got three datasets then plot them as a ternary system.

        #define the x/y system
        x = 0.5 - (concentrations[:,0]) * np.cos(np.pi/3) + (concentrations[:,1]/2)
        y = 0.866 - (concentrations[:,0]) * np.sin(np.pi/3) - ((concentrations[:,1] * (1/np.tan(np.pi/6))) / 2)

        #find coordinates of eutectic minimum
        idx = np.argmin(melt_points)
        xx = 0.5 - (concentrations[idx, 0]) * np.cos(np.pi/3) + (concentrations[idx,1]/2)
        yy = 0.866 - (concentrations[idx, 0]) * np.sin(np.pi/3) - ((concentrations[idx,1]*(1/np.tan(np.pi/6)))/2)

        plt.figure(figsize=(12,6))
        plt.subplot(121)

        build_ternary_axis(6, melt_points) # build axis for upper plot (melt)

        plt.axis('off') # get rid of the normal, boring axis.
        plt.axis('equal')
        plt.tricontourf(x, y, melt_points, 125, cmap='plasma')
        h = plt.colorbar(orientation='horizontal', label='Melting Point / °C', shrink=0.5)
        tick_locator = ticker.MaxNLocator(nbins=5)
        h.locator = tick_locator
        h.update_ticks()

        plt.plot([0,1,0.5,0],[0,0,np.sqrt(3)/2,0], color='black',linewidth='1') # outer triangle lines
        plt.scatter(xx, yy, facecolor='none', edgecolor='white') # put a circle at the eutectic point

        plt.subplot(122)

        build_ternary_axis(6, clear_points) # build axis for lower plot (clear)

        plt.axis('off') # get rid of the normal, boring axis.
        plt.axis('equal')
        plt.tricontourf(x, y, clear_points, 125, cmap='plasma')
        h = plt.colorbar(orientation='horizontal', label='Clearing Point / °C', shrink=0.5)
        tick_locator = ticker.MaxNLocator(nbins=5)
        h.locator = tick_locator
        h.update_ticks()

        plt.plot([0,1,0.5,0],[0,0,np.sqrt(3)/2,0], color='black',linewidth='1') # outer triangle lines
        plt.scatter(xx, yy, facecolor='none', edgecolor='white') # put a circle at the eutectic point
        plt.show()

    elif np.shape(concentrations)[1] == 4: # attempt a pyramid/heatmap plot

        print('Attempting pyramidal plot for ' + str(np.shape(concentrations)[1]) + ' components')

        fig = plt.figure(figsize=(10, 10), constrained_layout=True)
        ax1 = fig.add_subplot(121, projection='3d')

        plot_pyramidal_ax(ax1)
        label_points(ax1)
        plot_concs, plot_melts, plot_clear = get_random_subset(concentrations, melt_points, clear_points, 20)
        sct1 = plot_3d_tern(ax1, plot_concs, plot_melts, 0.05)
        cb1 = fig.colorbar(sct1, orientation='horizontal', label='Melting Point / °C', shrink=0.5)
        #cb1.set_alpha(1)
        #cb1.draw_all() # not needed in latest version of mpl?
        ax1.axis('off')


        ax2 = fig.add_subplot(122, projection='3d')
        plot_pyramidal_ax(ax2)
        label_points(ax2)
        sct2 = plot_3d_tern(ax2, plot_concs, plot_clear, 0.05)
        cb2 = fig.colorbar(sct2, orientation='horizontal', label='Clearing Point / °C', shrink=0.5)
        #cb2.set_alpha(1)
        #cb2.draw_all()
        ax2.axis('off')
        plt.show()

    else:
        print('Not plotting; ' + str(np.shape(concentrations)[1]) + ' components')

    return

# function to get user input
def get_input(prompt):
    return input(prompt).strip()

# function to parse input strings into lists of floats
def parse_input(input_str):
    return [float(x.strip()) for x in input_str.split(',') if x.strip()]

if __name__ == "__main__":
    '''
    if called in this way, prompt the user to give some inputs we can use
    '''
    melting_points_input = get_input("Please input melting points (°C), separated by commas: ")
    enthalpies_input = get_input("Please input enthalpies of fusion (kJ/mol), separated by commas: ")
    clearing_points_input = get_input("Please input clearing points (°C), separated by commas: ")

    melting_points = parse_input(melting_points_input)
    enthalpies = parse_input(enthalpies_input)
    clearing_points = parse_input(clearing_points_input)

    # ensure that all the lists have the same length
    if not (len(melting_points) == len(enthalpies) == len(clearing_points)):
        print("Error: All inputs must have the same number of components.")
        sys.exit(1)

    melting_points = np.array(melting_points)
    enthalpies = np.array(enthalpies)
    clearing_points = np.array(clearing_points)

    concentrations, melt_points, clear_points = make_phase_diagram(melting_points, enthalpies, clearing_points)

    print_composition(concentrations, melt_points, clear_points)

def plot_slice(concentrations, melt_points, component_index, component_value, alpha_value):
    slice_filter = (concentrations[:, component_index] >= component_value - 0.01) & (concentrations[:, component_index] <= component_value + 0.01)
    slice_concentrations = concentrations[slice_filter]
    slice_melt_points = melt_points[slice_filter]

    if len(slice_concentrations) < 3:
        print(f"Not enough points to plot for component slice at {component_value:.2f}")
        return None

    adjusted_concentrations = np.delete(slice_concentrations, component_index, axis=1)
    adjusted_concentrations /= np.sum(adjusted_concentrations, axis=1, keepdims=True)

    fig, ax = plt.subplots()
    ax.tricontourf(adjusted_concentrations[:, 0], adjusted_concentrations[:, 1], slice_melt_points, levels=15, cmap='viridis', alpha=alpha_value)
    plt.close(fig)  # Prevent showing the plot in the notebook/output

    return fig
    
def create_gif(concentrations, melt_points, filename='phase_diagram_slices.gif'):
    images = []
    for i in np.arange(0, 1.02, 0.02):  # Iterate over slices in 2% increments
        fig = plot_slice(concentrations, melt_points, 3, i, alpha_value=0.6)  # Adjust alpha_value as needed
        if fig is not None:
            fig.canvas.draw()  # Draw the figure to ensure it's fully rendered
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            images.append(image)
    
    if images:
        imageio.mimsave(filename, images, fps=30)  # Save the GIF
    else:
        print("No images to save.")  