import matplotlib.pyplot  as plt
import matplotlib.patches as patches
import numpy              as np
import scipy.signal       as signal
import scipy.stats        as stats
import bottleneck         as bn
import math

def smooth_z_tip(z_tip, method = "movmean"):
    '''
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
    dz_tip      = np.diff(z_tip)
    max_indx    = np.argmax(dz_tip)  # TODO: this isn't great. Should use peak find or get first zero crossing.
    if isinstance(max_indx, list): 
        max_indx = max_indx[0]
    start       = round(0.008*len(z_piezo))
    end         = max_indx# - round(start/2) # We want to be a little bit away from the peak
    
    # Find the first zero crossing before `end`
    while end - start > 10: # There has got to be a better way to do this.
        end -= 1
        if dz_tip[end] <= 0:
            break
    # tmp = dz_tip[start:end]
    # tmp = tmp[tmp <= 0]
    # if len(tmp) > 10:
    #     end = np.argmax(tmp)
    #     if isinstance(end, list): 
    #         end = end[-1]
    
    if end - start < 10:
        start = round(start/2)
    
    if end - start <= 10:
        start = 5
        end   = 20
    
    return (start, end)

def line_slope(z_piezo, z_tip, index = None):
    '''
    Algorithm for getting the slope of a force ramp. (The following
    diagram may not display properly in IDE tooltips.)
    
    \
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
     - `get_start_end()` takes the derivative of z_tip and finds the z_piezo location where that derivative is highest. This is the end point.
     - The start point is just 0.8% of the length of z_piezo (8 for ramps with 1024 samples; 4 for ramps with 512 samples)
     - The end point is reduced untile the first zero crossing (before the maximum) of the derivative of z_tip is found.
    
    This (start, end) value is used to fit to the linear region of the force ramp. If R^2 is greater than 0.9, then this slope is returned.
    Otherwise, decrease the end value, fit again, and check R^2. This is repeated until R^2 is greater than 0.9 or either of the following 
    condition is met:
     - There are less than 15 points between the start and end values.
     - The process has looped through 10 times without meeting either of the above criteria.
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
    # TM Deflection is called z_tip in the paper. Here I am using tm_defl to hold the entire 64x64 array of TM deflections.
    #tm_defl = data[:, 1, :]

    # Smooth the force ramp
    if smooth_func is not None:
        tm_defl = np.apply_along_axis(smooth_func, axis = 0, arr = tm_defl)

    # slope = -(tm_defl[e, :]-tm_defl[s, :]) / (z_piezo[e] - z_piezo[s])             # V / nm
    (slope, r2s, _, _) = line_slope(z_piezo, tm_defl)
    size               = int(math.sqrt(len(slope)))
    slope              = np.nan_to_num(slope)         # Get rid of NaN
    slope[slope <= 0]  = 0.000001                     # Get rid of zeros and negatives
    
    # Save?
    if savefile is not None and isinstance(savefile, str):
        #np.savetxt(savefile, slope, delimiter=",")
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
    mod_slope  = slope.copy()
    mod_slope[r2s < 0.9] = np.nan
    
    slice_s    = 0
    slice_e    = fixed_edge - 2
    left_slice = mod_slope[:, slice_s:slice_e]
    
    #fig, axs   = plt.subplots(1, 2, figsize=(8,4))

    # Plot the left part of the slopes
    # This average (highest bin in histogram) of left_slice is used to determine
    # the TM deflection sensitivity.
    #axs[0].imshow(left_slice)

    left_slice = left_slice.flatten()
    left_slice = left_slice[~np.isnan(left_slice)]  # Remove the points that had an R^2 < 0.9
    
    # Plot histogram of left_slice
    h            = np.histogram(left_slice.flatten(), bins=20)
    edges        = h[1]
    edges        = edges[1:]
    tm_defl_sens = 1/float(edges[h[0] == h[0].max()][-1])
            
    
    print(f"Sample = {sc_params['name']}")
    print(f"TM Defl. Sens. = {tm_defl_sens:.2f} nm/V")

    slope *= tm_defl_sens  # [V/nm]*[nm/V]=[1]


    comp            = 1/sc_params["afm_spring_constant"]*(slope**-1 - 1)     # Compliance
    comp[comp <= 0] = 0.000001

    if linearize:
        comp = comp**(1/3.0)
        
    return (comp, r2s)



# Move to new file?

def get_cantilever_params(params, cant_num):
    thick = params["thickness"]
    width = params["cantilevers"][cant_num]["width"]
    start = params["cantilevers"][cant_num]["start"]
    end   = params["cantilevers"][cant_num]["end"]
    igno  = params["cantilevers"][cant_num]["lin_ignore"]
    fixed = params["cantilevers"][cant_num]["fixed_edge"] - 1 # Parameters file indices start at 1, but in python they start at 0
    start = np.array(start) - 1
    end   = np.array(end)   - 1
    row   = (end[1] + start[1]) // 2 # Find center line of cantilever
    col_s = start[0] + igno
    col_e = end[0]
    
    return (thick, width, start, end, igno, fixed, start, end, row, col_s, col_e)

def comp_mat_inspector(comp_mat, z_piezo, tm_defl, params, fig_width = 12, r2s_mat = None):
    # Plot slopes images
    if 'fig' in locals():
        plt.close(fig)
    #fig, axs = plt.subplots(1, 2, figsize=(12,6))
    
    # fig, axs = plt.subplot_mosaic(
    #     [['left', 'upper right'],
    #      ['left', 'lower right']],
    #     figsize            = (fig_width, fig_width/2), 
    #     constrained_layout = True
    # )
    
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
        axs["D"].set_axis_off()
    

    size = int(math.sqrt(tm_defl.shape[1]))

    axs["A"].custom_info = {
        'z_piezo'   : z_piezo,
        'z_tip'     : tm_defl,
        'size'      : size,
        'ax_z_tip'  : axs["B"],
        'ax_dz_tip' : axs["C"]
    }
    
    axs["A"].pcolormesh(comp_mat, vmin=0, vmax=1)#, edgecolors='k', linewidth=0.1)
    axs["A"].invert_yaxis()
    axs["A"].set_title("Compliance map")
    
    if r2s_mat is not None:
        axs["D"].custom_info = axs["A"].custom_info
        
        axs["D"].pcolormesh(r2s_mat, vmin=0, vmax=1)#, edgecolors='k', linewidth=0.1)
        axs["D"].invert_yaxis()
        axs["D"].set_title("$R^2$ map")

    # Add lines over points to fit
    for cant_num in range(len(params["cantilevers"])):
        (thick, width, start, end, igno, fixed, start, end, row, col_s, col_e) = get_cantilever_params(params, cant_num)

        axs["A"].plot([col_s, col_e], [row+0.5, row+0.5], 'r')

        # Draw rectangle
        rect = patches.Rectangle((start[0], start[1]), end[0]-start[0], end[1]-start[1], linewidth=1, edgecolor='r', facecolor='none')

        axs["A"].add_patch(rect)

    plot_z_tip(0, 0, z_piezo, tm_defl, size, axs["B"], axs["C"])    
    cid = fig.canvas.mpl_connect('button_press_event', onclick_mat)
    
    return axs

def plot_z_tip(row, col, z_piezo, z_tip, size, ax1, ax2):
    '''
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

    #ax2.relim()
    #ax2.autoscale_view()
    plt.gcf().canvas.draw()
    #plt.gcf().canvas.flush_events()
    
    
    
def onclick_mat(event):
    '''
    '''
    try:
        z_piezo = event.inaxes.custom_info['z_piezo']
        z_tip   = event.inaxes.custom_info['z_tip']
        size    = event.inaxes.custom_info['size']
        ax1     = event.inaxes.custom_info['ax_z_tip']
        ax2     = event.inaxes.custom_info['ax_dz_tip']
    except:
        return

    col = math.floor(event.xdata)
    row = math.floor(event.ydata)
    col = 0  if col < 0  else col
    row = 0  if col < 0  else row
    col = (size - 1) if col > (size - 1) else col
    row = (size - 1) if row > (size - 1) else row
    
    # We use flipud on the matrix, so the row is wrong
    row = (size - 1) - row
    
    ax1.clear()
    ax2.clear()
    plot_z_tip(row, col, z_piezo, z_tip, size, ax1, ax2)