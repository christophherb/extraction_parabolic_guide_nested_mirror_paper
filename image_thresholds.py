import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import patches

def return_com_array(array: np.array):
    """short helper function to return the center of mass of an array in pixles

    Args:
        array (np.array): array of which to return the center of mass

    Returns:
        tuple: tuple(y, x) center of mass
    """
    height, width = array.shape[0], array.shape[1]
    X, Y = np.meshgrid(range(height), range(width))
    total = np.sum(array)
    y_mean, x_mean = np.sum(Y*array)/total, np.sum(X*array)/total
    return y_mean, x_mean

def return_thresholds(array, thresholds=None):
    """For a given list of thresholds and an intensity map,
    a list of (x, y) are returned such that the intensity contained by the centered rectangle with height 2*x and width 2*y
    surpasses the respective threshold.

    Args:
        array (np.array): array of intensities
        thresholds (list, optional): List of thresholds default is ten equidistant steps from 0.1 to 1. Defaults to None.

    Returns:
        list: list of tuple (x, y) for array-centered rectangle with heigth 2*x and width 2*y
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    threshold_xy = []
    shape = array.shape
    all = np.sum(array)
    midpoint = int(shape[0]/2), int(shape[1]/2)
    xmin = 0
    y_divided_x = round(shape[1]/shape[0]) # calc the y steps from x
    print(y_divided_x)
    for ind, threshold in enumerate(thresholds):
        for x in range(xmin, midpoint[0]):
            y = round(y_divided_x*x)
            try:
                data = array[midpoint[0]-x: midpoint[0]+x+1, midpoint[1]-y: midpoint[1]+y+1]
            except IndexError:
                break
            if np.sum(data)/all >= threshold:
                threshold_xy.append((x, y))
                xmin = x
                break
    return threshold_xy

def plot_thresholds(array, fig=None, ax=None, thresholds=None, extent=None, circle=None):
    """plots the calculated thresholds in a rectangular map

    Args:
        array (np.array): data array for wich the thresholds are calculated
        fig (plt.fig, optional): Figure in which to plot if none is given a new figure is created. Defaults to None.
        ax (plt.ax, optional): ax in which to plot if none is give a new one is created. Defaults to None.
        thresholds (list, optional): list of thresholds defaults to (0.1, 0.2, ..., 0.9). Defaults to None.
        circle (center, radius): if a center and radius are given

    Returns:
        (fig, ax): tuple of figure and ax
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    threshold_xy = return_thresholds(array, thresholds)
    midpoint = int(array.shape[0]/2), int(array.shape[1]/2)
    if extent==None:
        extent=(0, array.shape[0], 0, array.shape[1])
    height = extent[1]-extent[0]
    width = extent[3]-extent[2]

    def px2datapoints(px, py):
        x = (px/array.shape[0])*(height)+ extent[0]
        y = (py/array.shape[1])*(width)+ extent[2]
        return x, y

    if fig==None or ax==None:
        fig, ax = plt.subplots()
        im = ax.imshow(array[::-1], interpolation='none', extent=extent)

    for ind, (x, y) in enumerate(threshold_xy):
        lowerleft_x, lowerleft_y = px2datapoints(midpoint[0]-x, midpoint[1]-y)
        rect = patches.Rectangle((lowerleft_x, lowerleft_y), 2*x*height/array.shape[0], 2*y*width/array.shape[1], alpha=1, fc='none', edgecolor='red')
        ax.add_patch(rect)
        ax.text(lowerleft_x+0.2, lowerleft_y+0.2, round(thresholds[ind], 2), color='red')
    #fig.colorbar(im)
    return fig, ax


def return_thresholds_circle(array: np.array, xmid_px, ymid_px, extent: tuple=None, thresholds: list=None, normfunction=None) -> np.array:
    """finds and returns the radii of circles crossing certain thresholds

    Args:
        array (np.array): array of data to find thresholds for
        extent (tuple): extent of the array. Defaults to (0,xpix,0,ypix)
        thresholds (list): list of thresholds. Defaults to np.linspace(0.1, 0.9, 9)

    Returns:
        np.array: list(ind, threshold, r[px]) list for all thresholds
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    threshold_radii = []
    total = np.sum(array)
    radius = 0
    if normfunction is None:
        normfunction = lambda x, y, r: (x)**2+(y)**2 < r**2
    for ind, threshold in enumerate(thresholds):
        for r in range(radius, max([array.shape[0], array.shape[1]])):
            try:
                integrated_area = return_integrated_area(array, r, normfunction, xmid_px, ymid_px)
                if integrated_area/total > threshold:
                    threshold_radii.append((ind, threshold, r))
                    radius = r
                    break
            except IndexError:
                break
        else:
            print('no suc', threshold)
    return np.array(threshold_radii)

def plot_thresholds_circle(array: np.array, xmid_px=None, ymid_px=None, extent: tuple=None, thresholds: list=None, figax=None):
    if extent is None:
        extent = (-0.5, array.shape[0]-0.5, -0.5, array.shape[1]-0.5)
    if xmid_px is None:
        xmid_px = array.shape[0]//2
    if ymid_px is None:
        ymid_px = array.shape[1]//2
    midpointx, midpointy = extent[0] + (xmid_px/(array.shape[0]-1))*(extent[1]-extent[0]), extent[2] + (ymid_px/(array.shape[1]-1))*(extent[3]-extent[2])
    print(midpointx, midpointy)
    if figax is None:
        fig, ax = plt.subplots(1)
        im = ax.imshow(array[::-1], extent=extent, interpolation='none')
    else:
        fig, ax = figax
    # transform the coordinates to indices
    thresholds=return_thresholds_circle(array,  xmid_px=xmid_px, ymid_px=ymid_px, extent=extent, thresholds=thresholds)
    for ind, (i, threshold, r_px) in enumerate(thresholds):
        r_ext = r_px*(extent[1]-extent[0])/array.shape[0]
        print(r_ext, threshold)
        circ = patches.Circle((midpointx, midpointy), r_ext, fc='none', edgecolor='red')
        ax.add_patch(circ)
        ax.text(xmid_px/(array.shape[0]-1)+r_px/(array.shape[0]-1), ymid_px/(array.shape[1]-1),'$I_p$ = {:.2}'.format(threshold), color='red', transform=ax.transAxes)
    return fig, ax



def radius_vs_content(array: np.array, mid_x: float=None, mid_y: float=None, radii=None, normfunction=None):
    """takes an array and returns the radius of the cirle vs the incorporated content

    Args:
        array (np.array): array with data
        mid_x (float, optional): x-value of the midpoint. Defaults to None.
        mid_y (float, optional): y-value of the midpoint. Defaults to None.

    Returns:
        (np.array): [(rad, int)]
    """
    height = array.shape[0]
    width = array.shape[1]
    total = np.sum(array)
    rad_int = []
    if mid_x is None:
        mid_x = (height)/2
    if mid_y is None:
        mid_y = (width)/2
    if radii is None:
        radii = range(0, max([height, width]))
    if normfunction is None:
        normfunction = lambda x, y, r: (x)**2+(y)**2 <= r**2
    X, Y = np.meshgrid(range(height), range(width))
    for radius in radii:
        fraction = return_integrated_area(array, radius, normfunction, mid_x, mid_y)/total
        rad_int.append((radius, fraction))
    return rad_int

def return_integrated_area(array, r, function, xmid, ymid):
    height, width = array.shape
    X, Y = np.meshgrid(range(width), range(height))
    mask = norm((X-xmid), (Y-ymid),  r, function)
    return np.sum(array[mask])


def norm(X, Y, r, function):
    return function(X, Y, r)

a = np.ones((1000, 1000))
norm2 = lambda x, y, r: (x)**2+(y)**2 < r**2
print('this', return_integrated_area(a, 100, norm2, 500, 500))