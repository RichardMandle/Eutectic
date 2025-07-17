import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

def mpl_tri_plot(concentrations, melt_points, clear_points,cmap1='plasma',cmap2='viridis',figsize=6):
    '''
    Plot a ternary phase diagram using matplotlib
    
    concentrations  - np.array of concentrations (should be n-by-3 in size)
    melt_points     - np.array of melting points
    clear_points    - np.array of clearing points
    cmap1           - cmap to use for the left (melting point) plot
    cmap2           - cmap to use for the right (clearing point) plot 
    figsize         - size of the figure (will be 2*figsize x figsize)
    '''

    #define the x/y system
    x = 0.5 - (concentrations[:,0]) * np.cos(np.pi/3) + (concentrations[:,1]/2)
    y = 0.866 - (concentrations[:,0]) * np.sin(np.pi/3) - ((concentrations[:,1] * (1/np.tan(np.pi/6))) / 2)
    
    #find coordinates of eutectic minimum
    idx = np.argmin(melt_points)
    xx = 0.5 - (concentrations[idx, 0]) * np.cos(np.pi/3) + (concentrations[idx,1]/2)
    yy = 0.866 - (concentrations[idx, 0]) * np.sin(np.pi/3) - ((concentrations[idx,1]*(1/np.tan(np.pi/6)))/2)
    
    plt.figure(figsize=(2*figsize,figsize))
    plt.subplot(121)
    
    build_ternary_axis(6, melt_points) # build axis for upper plot (melt)
    
    plt.axis('off') # get rid of the normal, boring axis.
    plt.axis('equal')
    plt.tricontourf(x, y, melt_points, 125, cmap=cmap1)
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
    plt.tricontourf(x, y, clear_points, 125, cmap=cmap2)
    h = plt.colorbar(orientation='horizontal', label='Clearing Point / °C', shrink=0.5)
    tick_locator = ticker.MaxNLocator(nbins=5)
    h.locator = tick_locator
    h.update_ticks()
    
    plt.plot([0,1,0.5,0],[0,0,np.sqrt(3)/2,0], color='black',linewidth='1') # outer triangle lines
    plt.scatter(xx, yy, facecolor='none', edgecolor='white') # put a circle at the eutectic point
    plt.show()


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

    plt.text(0.5, -0.125, 'Component A / mol %', ha='center', va='center', fontsize=10)
    plt.text(0.1, 0.5, 'Component B / mol %', rotation=60, ha='center', va='center', fontsize=10)
    plt.text(0.9, 0.5, 'Component C / mol %', rotation=-60, ha='center', va='center', fontsize=10)


