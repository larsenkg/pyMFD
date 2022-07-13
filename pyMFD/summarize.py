from pyMFD.cantilever import get_cantilever_params

import matplotlib.pyplot  as plt
import matplotlib.patches as patches
import numpy              as np
import scipy.signal       as signal
import scipy.stats        as stats
import bottleneck         as bn
import math

def smooth_z_tip(z_tip, method = "movmean"):
    '''
    The raw force-deflection data from the AFM scan is frequently noisy. 
    This function performs either a moving average or a butterworth filter. 
    I found that doing two moving averages with a window size of 5 and then 
    11 works well. The data is first flipped to avoid a "lip" at the 
    beginning of the data.

    Parameter
    ---------
    z_tip : ndarray
        Numpy array containing the z_tip force-deflection data.
    method : {"movmean", "butter"}, optional, default: "movmean"
        Select smoothing method. Can be 'movmean' or 'butter'. 

    Returns
    -------
    z_tip_smooth : ndarray
        Filtered data. Same shape as `z_tip`.
    '''
    if method == "butter":
        b, a         = signal.butter(1, 0.1)
        z_tip_smooth = signal.filtfilt(b, a, z_tip)
    else: # movmean
        z_tip_smooth = bn.move_mean(np.flip(z_tip), window=5,  min_count=1, axis=0)
        z_tip_smooth = bn.move_mean(z_tip_smooth,   window=11, min_count=1, axis=0)
        z_tip_smooth = np.flip(z_tip_smooth)
        
    return z_tip_smooth

def get_start_end(z_piezo, z_tip):
    '''
    Get the start and end indices for the linear portion of z_piezo vs z_tip.

    Parameters
    ----------
    z_piezo : ndarray
        Z_piezo data as a numpy array.
    z_tip : ndarray
        Z_tip data as a numpy array.
    
    Returns
    -------
    start, end : int
        Start and end indices.
    '''
    dz_tip   = np.diff(z_tip)
    max_indx = np.argmax(dz_tip)  # TODO: this isn't great. Should use peak find or get first zero crossing.
    if isinstance(max_indx, list): 
        max_indx = max_indx[0]

    # The starting index should not be the very first, since the first few datapoints are noisy.
    # This is a silly way to do this, but we essentially discard the first 0.8% of the data.
    # This is 8 for samples with 1024 points and 4 for samples with 512 points.
    start = round(0.008*len(z_piezo))
    end   = max_indx   # This 'end' is a starting point for below
    
    # Find the first zero crossing before `end`
    while end - start > 10: # Repeat, but don't get to close to 'start'
        end -= 1
        if dz_tip[end] <= 0:
            break
    # tmp = dz_tip[start:end]
    # tmp = tmp[tmp <= 0]
    # if len(tmp) > 10:
    #     end = np.argmax(tmp)
    #     if isinstance(end, list): 
    #         end = end[-1]
    
    # If we started too close, then back off on 'start'
    # This will take 'start' in samples with 1024 points from 8 to 4.
    # This will take 'start' in samples with 512 points from 4 to 2.
    if end - start < 10:
        start = round(start/2)
    
    # As a last resort, just fix the start and end indices.
    if end - start <= 10:
        start = 5
        end   = 20
    
    return (start, end)

def line_slope(z_piezo, z_tip, index = None):
    '''
    Algorithm for getting the slope of a force ramp. Applied to all force 
    ramps in the force volume data (4096 for 64x64 scans). The slope is 
    used to find the compliance at each point in the map.

    (The following diagram may not display properly in IDE tooltips.)
    
    ::

        \\
         \                                  ^
          \                                 |
           \    __________________        z_tip
            \  /                            |
             \/                             v
        0123456789... <- z_piezo index
    
    
    Algorithm needs to find slope of the linear section from 0 to 5. 
    However, the data is rarely this nice. A robust algorithm is needed 
    to handle most cases.
    
    Get an initial (start, end) estimate using `get_start_end()`.
    
     - `get_start_end()` takes the derivative of z_tip and finds the 
       z_piezo location where that derivative is highest. This is the end 
       point.
     - The start point is just 0.8% of the length of z_piezo (8 for ramps 
       with 1024 samples; 4 for ramps with 512 samples)
     - The end point is reduced until the first zero crossing (before the 
       maximum) of the derivative of z_tip is found.
    
    This (start, end) value is used to fit to the linear region of the 
    force ramp. If R^2 is greater than 0.9, then this slope is returned.
    Otherwise, decrease the end value, fit again, and check R^2. This is 
    repeated until any of these conditions are met:

     - R^2 is greater than 0.9, or 
     - There are less than 15 points between the start and end values, or
     - The process has looped through 10 times without meeting either of 
       the above criteria.

    Parameters
    ----------
    z_piezo : ndarray
        Z_piezo data as a numpy array.
    z_tip : ndarray
        Z_tip data as a numpy array.
    index : int, optional
        If `index` is supplied, this function will not loop through all 
        force-deflection ramps in the FV data. It will only look at the 
        ramp where the index of z_tip is `index`. Useful for code that 
        selects only one force-ramp to plot.

    Returns
    -------
    slopes : ndarray
        Slopes for each force ramp in scan. Shape is (`size`,), where 
        `size` is the total number of force ramps in scan.
    r2s : ndarray
        R^2 array with same shape as `slopes`.
    s, e : int
        Start and end indices actually used to bracket region of interest.
    '''
    size   = z_tip.shape[1]       # 4096 for 64x64 scans; 1024 for 32x32 scans.
    slopes = np.zeros((size,)) 
    r2s    = np.zeros((size,))
    for n in range(size):
        if index is not None:
            n = index
            
        (s, e) = get_start_end(z_piezo, z_tip[:, n])
        r2     = 0
        count  = 0
        orig_e = e
        tries  = 10
        while r2 < 0.9:
            # Fit the region from `s` to `e`
            res = stats.linregress(z_piezo[s:e], z_tip[s:e, n])
            r2  = res.rvalue**2

            ## For next time
            old_e = e
            e    -= orig_e // tries
            
            # Make sure there are at least 20 samples
            if e - s < 15:
                e = old_e
                break
                
            count += 1
            if count >= tries:
                break
        
        slopes[n] = -res.slope  # Take the negative of slope
        r2s[n]    = r2
        
        # Only loop once if index is set to something
        if index is not None:
            break
        
    return (slopes, r2s, s, e)

def get_comp_mat(z_piezo, tm_defl, sc_params, linearize = True, savefile = None, smooth_func = smooth_z_tip, **kwargs):
    '''
    Get the compliance map. In other words, convert each force-deflection 
    ramp to a compliance value.

    Parameters
    ----------
    z_piezo : ndarray
        Piezo displacement data as a numpy array.
    tm_defl : ndarray
        Tapping mode deflection data as a numpy array.
    sc_params : dict
        Dictionary containg parameters loaded from JSON file with 
        `get_scan_params()`.
    linearize:  boolean, optional, default: True
        If true, will take the cube root of the compliance data. This 
        linearizes the data in displacement, since the compliance equation 
        depends on the position along the cantilever to the third power 
        (see Euler cantilever equation).
    savefile : str, optional
        If provided, the slopes will be saved to the file `savefile`.
    smooth_func : function, optional
        This function will be applied to `tm_defl` to smooth the 
        force-deflection data.

    Returns
    -------
    comp : ndarray
        Compliance matrix. Shape should be square, with the size of the 
        sides being the square root of the number of force ramps.
        E.g. shape is (64, 64).
    r2s : ndarray
        R^2 matrix. See `comp` for shape.
    '''
    # TM Deflection is called z_tip in the paper. Here I am using tm_defl 
    # to hold the entire 64x64 array of TM deflections.
    #tm_defl = data[:, 1, :]

    # Smooth the force ramp
    if smooth_func is not None:
        tm_defl = np.apply_along_axis(smooth_func, axis = 0, arr = tm_defl)

    (slope, r2s, _, _) = line_slope(z_piezo, tm_defl)
    size               = int(math.sqrt(len(slope)))
    slope              = np.nan_to_num(slope)         # Get rid of NaN
    slope[slope <= 0]  = 0.000001                     # Get rid of zeros and negatives
    
    # Save?
    if savefile is not None and isinstance(savefile, str):
        slope.tofile(savefile, sep=",", format="%2.8f")
        print(f"Saved to {savefile}.")
    
    slope  = slope.reshape((size, size))
    slope  = np.flipud(slope)
    r2s    = r2s.reshape((size, size))
    r2s    = np.flipud(r2s)
    
    # Get smallest fixed edge (closest to left edge of image)
    fixed_edge = size
    for cant in sc_params["cantilevers"]:
        if cant["fixed_edge"] < fixed_edge:
            fixed_edge = cant["fixed_edge"]

    # For determining the TM deflection sensitivity, ignore points with R^2 lower than 0.9
    mod_slope            = slope.copy()
    mod_slope[r2s < 0.9] = np.nan
    
    # Find TM deflection sensitivity by looking at the left part of the slope data.
    slice_s    = 0
    slice_e    = fixed_edge - 2
    left_slice = mod_slope[:, slice_s:slice_e]
    left_slice = left_slice.flatten()
    left_slice = left_slice[~np.isnan(left_slice)]  # Remove the points that had an R^2 < 0.9
    
    # Find highest bar in histogram of left_slice
    # Use that as the TM deflection sensitivity
    h            = np.histogram(left_slice.flatten(), bins=20)
    edges        = h[1]
    edges        = edges[1:]
    tm_defl_sens = 1/float(edges[h[0] == h[0].max()][-1])
    
    print(f"Sample = {sc_params['name']}")
    print(f"TM Defl. Sens. = {tm_defl_sens:.2f} nm/V")

    slope *= tm_defl_sens  # [V/nm]*[nm/V]=[1]

    # Calculate compliance using:
    #    comp = 1/k_afm*(1/slope - 1)
    comp            = 1/sc_params["afm_spring_constant"]*(slope**-1 - 1)     # Compliance
    comp[comp <= 0] = 0.000001 # Replace zeros with sufficiently small value (but not too small)

    if linearize:
        comp = comp**(1/3.0)
        
    return (comp, r2s)


def comp_mat_inspector(comp_mat, z_piezo, tm_defl, params, fig_width = 10, r2s_mat = None):
    '''
    Create the interactive compliance map inspector. This tool shows the 
    compliance map on the left, the selected force-deflection map in the 
    middle, and an R^2 map on the right. Click on any pixel in the 
    compliance map or the R^2 map to update the middle force-deflection 
    map.

    Parameters
    ----------
    comp_mat : ndarray
        Compliance matrix from `get_com_mat()`
    z_piezo : ndarray
        Piezo displacement data. Used for central plot.
    tm_defl : ndarray
        Tapping mode deflection data. Used for central plot.
    params : dict
        Dictionary of parameters. Load from JSON using `get_scan_params()`.
    fig_width : int, optional
        Width of matplotlib figure in inches.
    r2s_mat : ndarray, optional
        R^2 matrix to plot in third column. If not included, third column 
        is disabled.

    Returns
    -------
    axs : Axes
        Return matplotlibe axes used in figure.
    '''
    # Plot slopes images
    if 'fig' in locals():
        plt.close(fig)

    # Use matplotlib's mosaic to give nice names to axes.
    mosaic = """
        ABD
        ACD
        """
    
    fig = plt.figure(
        constrained_layout = True, 
        figsize            = (fig_width, fig_width/3)
    )
    axs = fig.subplot_mosaic(mosaic)
    
    if r2s_mat is None:
        # TODO: This turns the axis off, but the area is still there (just blank).
        axs["D"].set_axis_off()
    
    # Get number of pixels per row in map. Maps should always be square.
    size = int(math.sqrt(tm_defl.shape[1]))

    # Keep this info with the axis so that the interaction works.
    axs["A"].custom_info = {
        'z_piezo'   : z_piezo,
        'z_tip'     : tm_defl,
        'size'      : size,
        'ax_z_tip'  : axs["B"],
        'ax_dz_tip' : axs["C"]
    }
    
    # Plot compliance map
    axs["A"].pcolormesh(comp_mat, vmin=0, vmax=1)
    axs["A"].invert_yaxis()
    axs["A"].set_title("Compliance map")
    
    # Plot R^2 map
    if r2s_mat is not None:
        axs["D"].custom_info = axs["A"].custom_info
        
        axs["D"].pcolormesh(r2s_mat, vmin=0, vmax=1)
        axs["D"].invert_yaxis()
        axs["D"].set_title("$R^2$ map")

    # Add lines over points to fit
    for cant_num in range(len(params["cantilevers"])):
        (thick, width, igno, fixed, start, end, row, col_s, col_e) = get_cantilever_params(params, cant_num)

        axs["A"].plot([col_s, col_e], [row+0.5, row+0.5], 'r')

        # Draw rectangle
        rect = patches.Rectangle((start[0], start[1]), end[0]-start[0], end[1]-start[1], linewidth=1, edgecolor='r', facecolor='none')

        axs["A"].add_patch(rect)

    # Plot the z_tip data. Shows what was fit to get compliance.
    plot_z_tip(0, 0, z_piezo, tm_defl, size, axs["B"], axs["C"])    
    cid = fig.canvas.mpl_connect('button_press_event', onclick_mat)
    
    return axs

def plot_z_tip(row, col, z_piezo, z_tip, size, ax1, ax2):
    '''
    Plot the z_tip data. Each pixel in the compliance map comes from 
    fitting to z_tip.

    Parameters
    ----------
    row : int
        Row from compliance map. Used along with `col` to identify specific 
        pixel.
    col : int
        Column from compliance map. Used along with `row` to identify s
        pecific pixel.
    z_piezo : ndarray
        Piezo displacement data.
    z_tip : ndarray
        AFM tip displacement data.
    size : int
        Number of columns per compliance map.
    ax1, ax2 : Axes
        Two axes on which to plot. `ax1` is used for the z_tip data and 
        `ax2` is used of its derivative.
    '''
    index  = row*size + col
    ax1.plot(z_piezo, z_tip[:, index]*1000)
    ax1.set_xlabel("$Z_{piezo}$ (nm)")
    ax1.set_ylabel("$Z_{tip}$ (mV)")
    
    (_, _, s, e) = line_slope(z_piezo, z_tip, index = index)
    x            = np.array(range(0, z_tip.shape[1]))
    res          = stats.linregress(z_piezo[s:e], z_tip[s:e, index]*1000)
    print(f"({col}, {row}) : {index}; slope = {res.slope:.3f}; tm sens. = {-1/res.slope*1000:.1f} nm/V; $r^2$ = {res.rvalue**2:.2f}")
    ax1.set_title(f"({col}, {row}) : {index}; slope = {res.slope:.3f}; tm sens. = {-1/res.slope*1000:.1f} nm/V; $r^2$ = {res.rvalue**2:.2f}")
    ax1.plot(z_piezo[s:e], res.intercept + res.slope*z_piezo[s:e])
    
    # Plot derivative
    ax2.plot(z_piezo[1:], np.diff(z_tip[:, row*size + col]*1000))
    ax2.set_xlabel("$Z_{piezo}$ (nm)")
    ax2.set_ylabel("$dZ_{tip}$ (mV)")
    
    ax2.axvline(x = z_piezo[s], c='c')
    ax2.axvline(x = z_piezo[e], c='m')

    plt.gcf().canvas.draw()
    
    
    
def onclick_mat(event):
    '''
    Click event hander. Used to allow for inspection of the compliance map.

    Parameters
    ----------
    event : matplotlib.backend_bases.Event
        Event fired when mouse clicked on compliance map.
    '''
    try:
        z_piezo = event.inaxes.custom_info['z_piezo']
        z_tip   = event.inaxes.custom_info['z_tip']
        size    = event.inaxes.custom_info['size']
        ax1     = event.inaxes.custom_info['ax_z_tip']
        ax2     = event.inaxes.custom_info['ax_dz_tip']
    except:
        return

    # Find the column and row of the pixel that was clicked
    col = math.floor(event.xdata)
    row = math.floor(event.ydata)
    col = 0  if col < 0  else col
    row = 0  if col < 0  else row
    col = (size - 1) if col > (size - 1) else col
    row = (size - 1) if row > (size - 1) else row
    
    # We use flipud on the matrix, so the row is wrong
    row = (size - 1) - row
    
    # Replot
    ax1.clear()
    ax2.clear()
    plot_z_tip(row, col, z_piezo, z_tip, size, ax1, ax2)