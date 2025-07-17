import numpy as np
from itertools import combinations
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mayavi import mlab

'''
get_verts                               - returns a nd.array of the vertices needed to make a tetrahedral pyramid (no arguments)
plot_pyramidal_ax                       - plots the tetrahedral axis using mlab
plot_pyramidal_ax_mpl                   - plots the tetrahedral axis using matplotlib
label_pyramidal_axes                    - add tick labels to the pyramidal axis using mlab
mvi_quat_scatter                        - do a mlab scatter plot of the quaternary phase diagram
mvi_quat_isosurf                        - do a mlab isosurface plot of the quaternary phase diagram
label_points                            - labels the pyramidal axis
get_cartesian_array_from_barycentric    - gets a cartesian array from the barycentre of the pyramidal plot
plot_3d_tern                            - plots a 3D scatter plot within a pyramidal axis
get_random_subset                       - generates a random subset of compositions which includes the eutectic

rough code for quaternary phase diagram plotting with mayavi...
'''

def get_verts():
    '''
    This just generates and returns the verteces that define a pyramid.
    '''
    return np.array([
        [0, 0, 0],                           # Vertex 1
        [1, 0, 0],                           # Vertex 2
        [0.5, np.sqrt(3) / 2, 0],            # Vertex 3
        [0.5, 0.28867513, 0.81649658]        # Vertex 4
    ])

def plot_pyramidal_ax():
    '''
    Just draws the outline of the pyrmid using plot3d
    '''
    verts = get_verts()
    for v1, v2 in combinations(verts, 2):  # All edges
        line = np.array([v1, v2])
        mlab.plot3d(line[:, 0], line[:, 1], line[:, 2],
                    color=(0.5, 0.5, 0.5), tube_radius=0.005)

def label_pyramidal_axes():
    '''
    Handles tick labels for the pyramid
    '''
    ticks = np.linspace(0, 1, 5) # to change number of ticks, edit here (could expose this)
    verts = get_verts()
    for t in ticks:
        mlab.text3d(t, 0, 0, f'{t:.2f}', scale=0.05, color=(0, 0, 0))
        mlab.text3d(1 - t / 2, t * np.sqrt(3) / 2, 0, f'{t:.2f}', scale=0.05, color=(0, 0, 0))
        mlab.text3d(t / 2, t * np.sqrt(3) / 2, 0, f'{t:.2f}', scale=0.05, color=(0, 0, 0))
        mlab.text3d(0.5, 0.28867513, 0.81649658 * t, f'{t:.2f}', scale=0.05, color=(0, 0, 0))

def mvi_quat_scatter(c, m, i):
    '''
    scatter plot for quaternary phase diagram.

    input:
        c: np.ndarray, concentrations (X by 4 array).
        m: np.ndarray, scalar values corresponding to the concentrations.
    '''
    verts = get_verts()  # grab verteces
    x, y, z = np.dot(c, verts[:, 0]), np.dot(c, verts[:, 1]), np.dot(c, verts[:, 2])

    mlab.figure(size=(800, 600), bgcolor=(1, 1, 1))

    plot_pyramidal_ax()
    label_pyramidal_axes()

    mlab.points3d(x, y, z, m,
                  mode='sphere',
                  colormap='viridis',
                  scale_mode='none',
                  scale_factor=0.02)

    mlab.colorbar(title='Melting Point', orientation='vertical')

    mlab.show()

def mvi_quat_isosurf(c, m, i):
    '''
    quaternary isosurface using mayavi.

    input:
        c: np.ndarray, concentrations (X by 4 array).
        m: np.ndarray, scalar values corresponding to the concentrations.
    '''
    verts = get_verts()  
    x, y, z = np.dot(c, verts[:, 0]), np.dot(c, verts[:, 1]), np.dot(c, verts[:, 2])

    grid_x, grid_y, grid_z = np.mgrid[0:1:100j, 0:np.sqrt(3) / 2:100j, 0:1:100j]
    grid_values = griddata((x, y, z), m, (grid_x, grid_y, grid_z), method='linear', fill_value=np.nan)

    mlab.figure(size=(800, 600), bgcolor=(1, 1, 1))

    plot_pyramidal_ax()
    label_pyramidal_axes()

    scalar_field = mlab.pipeline.scalar_field(grid_x, grid_y, grid_z, grid_values)
    mlab.pipeline.iso_surface(scalar_field,
                              contours=[np.nanmean(m)],
                              opacity=0.5,
                              colormap='viridis')

    mlab.colorbar(title='Melting Point', orientation='vertical')

    mlab.show()

def mpl_quat_plot(c,m,i):
    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    ax1 = fig.add_subplot(121, projection='3d')
    
    plot_pyramidal_ax_mpl(ax1)
    label_points(ax1)
    plot_concs, plot_melts, plot_clear = get_random_subset(c, m, i, 20)
    sct1 = plot_3d_tern(ax1, plot_concs, plot_melts, 0.05)
    cb1 = fig.colorbar(sct1, orientation='horizontal', label='Melting Point / °C', shrink=0.5)
    #cb1.set_alpha(1)
    #cb1.draw_all() # not needed in latest version of mpl?
    ax1.axis('off')
    
    ax2 = fig.add_subplot(122, projection='3d')
    plot_pyramidal_ax_mpl(ax2)
    label_points(ax2)
    sct2 = plot_3d_tern(ax2, plot_concs, plot_clear, 0.05)
    cb2 = fig.colorbar(sct2, orientation='horizontal', label='Clearing Point / °C', shrink=0.5)
    #cb2.set_alpha(1)
    #cb2.draw_all()
    ax2.axis('off')
    plt.show()

def get_random_subset(concentrations, melt_points, clear_points, n):
    """
    Obtain a random subset of points and find the eutectic values.

    The function selects a random subset of points from the given
    concentration array (concentrations) and corresponding melting points array (melt_points),
    based on the specified subsampling interval (n). It also identifies the
    lowest value in the subset of melting points and its corresponding
    concentration values.

    random choice will always contain the eutectic (will it?)

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
    
def plot_pyramidal_ax_mpl(axis):
    """
    Plot the tetrahedral outline in 3D.

    The function plots the tetrahedral outline defined by the vertices
    using lines connecting the vertices.

    No inputs, no returns
    """

    verts = get_verts()
    
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
    
    verts = get_verts()
    
    t = np.transpose(np.array(verts))
    t_array = np.array([t.dot(x) for x in b])
    return t_array
