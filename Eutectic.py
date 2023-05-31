import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import ticker
from matplotlib.patches import Circle, PathPatch
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

import string
import os

'''
Simple python script for calculating phase diagrams and eutectic points.

Functions:
SaveData                                - saves data to .csv file (warning, files can be huge)
BuildTernaryAxis                        - called during plotting; makes the ternary axis system used

plot_pyramidal_ax                       - generates a triangular pyramidal axis for 4-component plotting
label_points                            - labels the pyramidal axis
get_cartesian_array_from_barycentric    - gets a cartesian array from the barycentre of the pyramidal plot
plot_3d_tern                            - plots a 3D scatter plot within a pyramidal axis
get_random_subset                       - generates a random subset of compositions which includes the eutectic

PrintComposition                        - prints the composition of the eutectic point (default) or ...
DoPhaseDiagram                          - Phase diagram engine - calculates the melting point at each concentration
MakePlot                                - Plotting engine - makes plots for 2- or 3- component systems.

'''

def SaveData(Concs,Melt,Clear,Filename='EutecticData.csv'):
    '''
    Save data to CSV file
    Headers tell you whats what.
    '''
    if os.path.exists(Filename): # check if Filename exists
        # If it does, add a suffix to the Filename to make it unique
        suffix = 1
        while True:
            new_filename = "{} ({}){}".format(*os.path.splitext(Filename) + (suffix,))
            if not os.path.exists(new_filename):
                Filename = new_filename
                break
            suffix += 1
            
    Headers = ["conc {}".format(i) for i in range(1,Concs.shape[1]+1)]
    Headers.extend(["Melt", "Clear"])

    data = np.column_stack((Concs,Melt,Clear))     # Stack the arrays horizontally

    np.savetxt(Filename, data, delimiter=",", header=",".join(Headers), fmt="%s")     # Save the data to a CSV file
    
    print('saved data to: ' + Filename)
    
    return
    
    
def BuildTernaryAxis(Divisions,ZValues):
    # Builds the triangular axis system used for ternary plots
    
    d1 = np.cos(np.pi/3)
    d2 = np.sin(np.pi/3)

    l = np.linspace(0,1,Divisions)
    zmax = max(ZValues)
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

    Parameters:
    None

    Returns:
    None
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

    Parameters:
    None

    Returns:
    None
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


def plot_3d_tern(axis,concs, melting_points,alpha_array):
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
    Scatter = axis.scatter(cartesian_points[:-1, 0], cartesian_points[:-1, 1], cartesian_points[:-1, 2],
               c=melting_points[:-1], cmap='magma', s=100, alpha=alpha_array)
    
    #do the eutectic point slightly differently.
    axis.scatter(cartesian_points[-1, 0], cartesian_points[-1, 1], cartesian_points[-1, 2],
               c=melting_points[-1], cmap='magma', s=100, alpha=1)
    
    return(Scatter)

    
def get_random_subset(Concs, Melt, Clear, n):
    """
    Obtain a random subset of points and find the eutectic values.

    The function selects a random subset of points from the given
    concentration array (Concs) and corresponding melting points array (Melt),
    based on the specified subsampling interval (n). It also identifies the
    lowest value in the subset of melting points and its corresponding
    concentration values.
    
    random choice will always contain the eutectic

    Parameters:
    - Concs (ndarray): Array of concentration values.
    - Melt (ndarray): Array of melting points.
    - n (int): Subsampling interval.

    Returns:
    - random_concs (ndarray): Random subset of concentration values.
    - random_melts (ndarray): Random subset of melting points.

    """
    
    random_indices = np.append(np.random.choice(range(Concs.shape[0]), size=Concs.shape[0] // n, replace=False),
                               np.argmin(Melt))
    random_concs = np.squeeze(Concs[random_indices,:])
    random_melts = np.squeeze(Melt[random_indices])
    random_clear = np.squeeze(Clear[random_indices])

    return random_concs, random_melts, random_clear


def PrintComposition(Concs,Melts,Clear):
    '''
    Simply prints the eutectic composition, i.e. the lowest melt, from an array of concentrations and corresponding melting/clearing points
    
    Inputs:
    Concs - an array of concentrations for each mixtre
    Melts - a list of melting points for each mixture (in K)
    Clear - a list of clearing points for each mixture (in K)
    
    returns:
    some text
    '''

    print ('eutectic composition:')
    for n in range (0,np.shape(Concs)[1]):
        print (str(np.double(np.round(Concs[np.where(Melts==np.min(Melts)),n],decimals=3))) + ' mol% ' + chr(ord('A')+n))
    
    print ('Melting point: ' + str(np.round(np.min(Melts),decimals=2)) + ' °C')
    print ('Clearing point: ' + str(np.double(np.round(Clear[np.where(Melts==np.min(Melts))],decimals=2))) + ' °C')
    print ('enantiotropic range = ' + str(np.round(np.double(Clear[np.where(Melts==np.min(Melts))]-np.min(Melts)),decimals=2)) + ' °C')
    
    return
    

def DoPhaseDiagram(MeltPoints,ClearingPoints,EnthalpFusion,SearchSize=250000,Plotting=True):
    '''
    Takes a list of melting points, enthalpies of fusion, molecular weights
    searches a number of random concentrations (set by SearchSize) and calculates
    the Melting point and Clearing point of each.
    
    Works well for non-LCs & nematics, less well for smectics, untested for everything else.
    
    Basic idea is taken from https://doi.org/10.1039/C39740000098 but implemented differently;
    rather than seek the absolute eutectic we just sample a huge number of different concentrations.
    
    Inputs:
    MeltPonts           - the melting point of each component in Kelvin
    ClearingPoints      - the clearing point of each component in Kelvin
    EnthalpFusion       - the enthalpy of fusion of each component (i.e. melt) in kJ mol-1
    SearchSize (250000) - the number of concentrations to search
    Plotting (True)     - wether or not to make a plot; this does slow it down a bit.
    
    Returns:
    Concs - the concentrations used; a NxM array, where N is the number of components and M is the SearchSize
    Melt  - the melting point of each different concentration
    Clear - the clearing point of each different concentration
    '''

    Temp = np.random.rand(SearchSize,np.size(MeltPoints))
    Concs = Temp/np.transpose(np.tile(np.sum(Temp,1),(np.size(MeltPoints),1)))

    M = np.double(np.tile(MeltPoints,(SearchSize,1)))
    E = np.double(np.tile(EnthalpFusion,(SearchSize,1)))
    C = np.double(np.tile(ClearingPoints,(SearchSize,1)))

    Melt = (np.max((E*1000)/(((E*1000)/(M+273.15)-(8.314*np.log(Concs)))),1))-273.15
    Clear = np.sum(C*Concs,1)
    
    if Plotting:
        MakePlot(Concs,Melt,Clear)

    return(Concs,Melt,Clear)
    
    
def MakePlot(Concs,Melt,Clear):
    '''
    Plotting engine. Will attempt to make a plot; if we have 2x components it'll be a line graph, if 3 a ternary one.
    
    Inputs:
    Concs - array of concentrations of each component in a given mixture
    Melt  - array of melting points of each mixture
    Clear - array of clearing points of each mixture
    
    Returns:
    A nice plot (if you've got a 2- or 3- component mix!) or a warning that you have too many things to plot.
    '''
    
    if np.shape(Concs)[1]==2:
        fig = plt.figure()
        MeltLine = plt.scatter(Concs[:,1],Melt,label='Melting Point / °C')
        ClearLine = plt.scatter(Concs[:,1],Clear,label='Clearing Point / °C')
        plt.xlabel('Concentration of B / mol %')
        plt.ylabel('Temperature / °C')
        plt.legend(loc='upper right')
        plt.xlim([0,1])

        # stick a line on the eutectic composition
        X = np.ravel([Concs[np.where(Melt==np.min(Melt)),1],Concs[np.where(Melt==np.min(Melt)),1]])
        Y = np.ravel([np.min(Melt),Clear[np.where(Melt==np.min(Melt))]])
        plt.plot(X,Y,color='black',linewidth='0.5',linestyle='dashed')
        plt.show()
                
    # if we've got three datasets then plot them as a ternary system.
    if np.shape(Concs)[1]==3:

        #define the x/y system
        x = 0.5-(Concs[:,0])*np.cos(np.pi/3)+(Concs[:,1]/2)
        y = 0.866-(Concs[:,0])*np.sin(np.pi/3)-((Concs[:,1]*(1/np.tan(np.pi/6)))/2)

        #find coordinates of eutectic minimum
        xx = 0.5-(Concs[np.where(Melt==np.min(Melt)),0])*np.cos(np.pi/3)+(Concs[np.where(Melt==np.min(Melt)),1]/2)
        yy = 0.866-(Concs[np.where(Melt==np.min(Melt)),0])*np.sin(np.pi/3)-((Concs[np.where(Melt==np.min(Melt)),1]*(1/np.tan(np.pi/6)))/2)

        plt.figure(figsize=(12,6))
        plt.subplot(121)

        BuildTernaryAxis(6,Melt) # build axis for upper plot (melt)

        plt.axis('off') # get rid of the normal, boring axis.
        plt.axis('equal')
        plt.tricontourf(x,y,Melt,125,cmap='plasma')
        h = plt.colorbar(orientation='horizontal', label='Melting Point / °C',shrink=0.5)
        tick_locator = ticker.MaxNLocator(nbins=5)
        h.locator = tick_locator
        h.update_ticks()

        plt.plot([0,1,0.5,0],[0,0,np.sqrt(3)/2,0], color='black',linewidth='1') # outer triangle lines
        plt.scatter(xx,yy,facecolor='none',edgecolor='white') # put a circle at the eutectic point

        plt.subplot(122)

        BuildTernaryAxis(6,Clear) # build axis for lower plot (clear)

        plt.axis('off') # get rid of the normal, boring axis.
        plt.axis('equal')
        plt.tricontourf(x,y,Clear,125,cmap='plasma')
        h = plt.colorbar(orientation='horizontal',label='Clearing Point / °C',shrink=0.5)
        tick_locator = ticker.MaxNLocator(nbins=5)
        h.locator = tick_locator
        h.update_ticks()

        plt.plot([0,1,0.5,0],[0,0,np.sqrt(3)/2,0], color='black',linewidth='1') # outer triangle lines
        plt.scatter(xx,yy,facecolor='none',edgecolor='white') # put a circle at the eutectic point
        plt.show()
        
    if np.shape(Concs)[1] ==4:
        # attempt a pyramid/heatmap plot
        print('Attempting pyramidal plot for ' + str(np.shape(Concs)[1]) + ' components')
        
        fig = plt.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(121, projection='3d')

        plot_pyramidal_ax(ax1)
        label_points(ax1)
        plot_concs,plot_melts,plot_clear = get_random_subset(Concs,Melt,Clear,20)
        sct1 = plot_3d_tern(ax1,plot_concs, plot_melts,0.05)
        cb1 = fig.colorbar(sct1,orientation='horizontal',label='Melting Point / °C',shrink=0.5)
        cb1.set_alpha(1)
        cb1.draw_all()
        ax1.axis('off')


        ax2 = fig.add_subplot(122, projection='3d')
        plot_pyramidal_ax(ax2)
        label_points(ax2)
        sct2 = plot_3d_tern(ax2,plot_concs, plot_clear,0.05)
        cb2 = fig.colorbar(sct2,orientation='horizontal',label='Clearing Point / °C',shrink=0.5)
        cb2.set_alpha(1)
        cb2.draw_all()
        ax2.axis('off')
        plt.show()
        
        
    # if we've got 5 just print a warning
    if np.shape(Concs)[1]>=5:
        print('Not plotting; ' + str(np.shape(Concs)[1]) + ' components')
        
            
    return
    
    