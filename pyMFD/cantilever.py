import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as stats

def get_cantilever_params(params, cant_num):
    '''
    Get the important parameters from the parameter dictionary (loaded from JSON) for a specific cantilever.

    Parameters
    ----------
    params: dict
        Dictionary of parameters. Load from JSON using `get_scan_params()`. Pass in only parameter for single sample.
    cant_num: int
        Cantilever number for which to get params.
    '''
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
    
    return (thick, width, igno, fixed, start, end, row, col_s, col_e)

def get_cantilever_pos(pixel_size, size):
    '''
    Returns the pixel locations in meters across a row.

    Parameters
    ----------
    pixel_size: float
        Size of pixel in meters.
    size: int
        Size of scan (number of pixels in scan).
    '''
    return pixel_size*np.arange(0, size, 1)

def get_compliance_row(comp_mat, row, rows_to_avg = 1):
    '''
    Return a full row of the compliance map. If rows_to_avg is greater than 1, then rows above and below `row` will be averaged.

    Parameters
    ----------
    comp_mat: ndarray
        Compliance matrix
    row: int
        Row of scan to extract
    rows_to_avg: int, optional
        Total number of rows to average. Will always be symmetric, rounded up. Passing in 2 or 3 is equivalent
    '''
    one_sided = rows_to_avg // 2
    comp_row  = comp_mat[(row - one_sided):(row + one_sided + 1), :]   
    comp_row  = np.nanmean(comp_row, axis = 0)
    
    return comp_row

def fit_compliance_linear(position, compliance):
    '''
    Fit the linearized position vs compliance graph.
    
    1/k       = (4/(E*w*t^3))*(L-c)^3
    1/k^(1/3) = (4/(E*w*t^3))^(1/3)*(L-c)
    1/k^(1/3) = a*(L-c)
    1/k^(1/3) = a*L-a*c
    y   = m*x+b
    
    The slope is proportional to E. To get the fixed end offset, divide the intercept by the slope (and take negative).

    Parameters
    ----------
    position: ndarray
        Vector of positions (in meters).
    compliance: ndarray
        Vector of linearized compliance.
    '''
    res = stats.linregress(position, compliance)
    return (res.slope, res.intercept)

def fit_fun(L, a, c):
    '''
    Function to fit with scipy.optimize.curve_fit.

    Parameters
    ----------
    L: float
        Position along cantilever
    a: float
        Combination of width, thickness, and modulus.
    c: float
        Offset from L.
    '''
    return 1/a*(L - c)**3

def standardize_and_fit(positions, compliances, width, thickness, func = fit_fun):
    '''
    Standardized and then fit the non-linearized (i.e. original) compliance data.
    Parameters
    ----------
    position: ndarray
        Vector of positions (in meters).
    compliance: ndarray
        Vector of linearized compliance.
    width: float
        Width of cantilever
    thick: float
        Thickness of cantilever
    func: function

    '''

    # Remove mean and set standard deviation to 1
    mu      = np.mean(positions)
    sigma   = np.std(positions)
    pos_std = (positions - mu) / sigma
    
    popt, _  = curve_fit(fit_fun, pos_std, compliances, [150, -1.5])
    Y        = 4*popt[0]*sigma**3/(width*thickness**3)
    pos_off  = popt[1]*sigma + mu
    return (Y, pos_off, popt[0]*sigma**3)
    
def calc_modulus_offset(slope, intercept, width, thickness):
    '''
    Calculate the modulus and fixed end offset.

    Parameters
    ----------
    slope: float
        Slope returned from `fit_compliance_linear()`
    intercept: float
        Intercept returned from `fit_compliance_linear()`
    width: float
        Width of cantilever
    thick: float
        Thickness of cantilever
    '''
    E = 4/(slope**3 * width * thickness**3)
    c = -intercept / slope
    
    return (E, c)

def offset_to_col_coord(offset, col_s, pixel_size):
    '''
    TODO: Think of a better name for this function.
    TODO: Remove if unused
    
    Takes the offset calculated from fitting the the compliance row and return where that offset is in compliance map space.

    Parameters
    ----------
    offset: float
        Offset calculated by fitting. This is how far the estimated fixed end is from the origin selected for `position` array.
    col_s: 
    pixel_size: float
        Size of single pixel in meters.
    '''
    return round(col_s + offset / pixel_size)