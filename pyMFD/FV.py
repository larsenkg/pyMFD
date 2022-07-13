from pyMFD.nanoscope   import get_fv_data, get_params
from pyMFD.scan_params import get_scan_params
from pyMFD.summarize   import get_comp_mat

class FV:
    '''
    This class represents a single force-volume scan. It contains the 
    relevant scan parameters and force-volume data.

    Attributes
    ----------
    fv_filename : str
        String pointing to the force-volume scan file.
    sc_params_filename : str
        String pointing to the force-volume scan parameters file.
    fv_params : dict
        Parameters extracted from the force-volume scan file header.
    sc_params : dict
        Scan parameters extracted from scan parameters JSON file.
    z_piezo : ndarray
        Displacement of AFM piezo. Has size `ramp_len` (from parameters 
        `samples_per_ramp`).
    tm_defl : ndarray
        Tapping mode deflection. Has shape (`ramp_len`, 2, `num_curves`). 
        The 2 comes from having both an extension and retraction and 
        `num_curves` is the number of force-ramps in a scan, e.g. 4096.
    pixel_size : float
        Size of single pixel in force-volume map, in meters.
    _fv_params_func : function
        Function that returns the `fv_params` dictionary.
    _fv_data_func : function
        Function that returns the force-volume data.
    _sc_params_func : function
        Function that returns the `sv_params` dictionary.

    Methods
    -------
    get_extend(): 
        Get extenstion ramp data
    get_retract():
        Get retraction ramp data
    get_pixel_size([scan_size, scan_points]): 
        Get pixel size
    summarize([which_dir, summary_func]):
        Summarize the ramp data (i.e. extract compliance)
    '''

    def __init__(
        self, 
        fv_filename, 
        fv_params_func     = get_params, 
        fv_data_func       = get_fv_data, 
        sc_params_filename = None, 
        sc_params_func     = get_scan_params,
        fv_params_kwargs   = {},
        fv_data_kwargs     = {},
        sc_params_kwargs   = {},
        ):
        '''
        Create a new FV class.
        
        Parameters
        ----------
        fv_filename : str
        fv_params_func : function, optional
            Function that takes a string to the FV file and returns the 
            required parameters in a dictionary. See `pyMFD.nanoscope` for 
            nformation on the required parameters.
        fv_data_func : function, optional
            Function that takes a string to the FV file and parameters 
            dictionary. Must return a tuple where the first element is a 
            1-D np.ndarray containing the z_piezo ramp deflection series. 
            The second argument is a np.ndarray with shape (X, Y, Z) 
            containing the force-volume tip deflection data. X should be 
            the size of the 1-D z_piezo ramp, Y should be 1 or 2 
            (depending on if only trace, or trace and retrace are 
            included), and Z should be the squared value of the size of the 
            FV scan. E.g. Z=4096 for a 64x64 "pixel" scan.
        sc_params_filename : str, optional
            String containing the path to the scan parameters filename.
        sc_params_func : function, optional
            Function that takes sc_params_filename and returns a dictionary 
            containing the required scan parameters. See pyMFD.scan_params 
            for information on the required parameters.
        '''
        self.fv_filename             = fv_filename
        self.sc_params_filename      = self.fv_filename + ".json" if sc_params_filename is None else sc_params_filename
        self._fv_params_func         = fv_params_func
        self._fv_data_func           = fv_data_func
        self._sc_params_func         = sc_params_func
        self.fv_params               = self._fv_params_func(self.fv_filename, **fv_params_kwargs)
        (self.z_piezo, self.tm_defl) = self._fv_data_func(self.fv_filename, self.fv_params, **fv_data_kwargs)
        self.sc_params               = self._sc_params_func(self.sc_params_filename, **sc_params_kwargs)
        self.pixel_size              = self.get_pixel_size()


    def get_pixel_size(self, scan_size=None, scan_points=None):
        '''
        Calculate the size of a single pixel in the force-volume data. 
        Should be in units of meters.

        Parameters
        ----------
        scan_size : float, optional
            The total size of the force-volume scan (in meters).
        scan_points : int, optional
            The number of force-deflection ramps in each line of the scan.

        Returns
        -------
        float
            Size of pixel (in meters).
        '''
        if scan_size is None:
            scan_size   = self.fv_params["scan_size"]

        if scan_points is None:
            scan_points = self.fv_params["ramps_per_line"]

        return scan_size / scan_points

    def get_extend(self):
        '''
        Return the force-volume data recorded during the extension of the 
        AFM cantilever.

        Returns
        -------
        ndarray
            The extension curves of the tapping mode deflection data.
            Shape is (ramp_length, 1, num_curves), e.g. (1024, 1, 4096)
        '''
        return self.tm_defl[:, 0, :]

    def get_retract(self):
        '''
        Return the force-volume data recorded during the retraction of the 
        AFM cantilever.

        Returns
        -------
        ndarray
            The retraction curves of the tapping mode deflection data.
            Shape is (ramp_length, 1, num_curves), e.g. (1024, 1, 4096)
        '''
        return self.tm_defl[:, 1, :]

    def summarize(self, which_dir='retrace', summary_func=get_comp_mat, **kwargs):
        '''
        Create a 2D representation of the force-volume data.
        
        Parameters
        ----------
        which_dir : str {'trace', 'extend', 'retrace', 'retract'}, int {0 for trace, 1 for retrace}
            Select whether the trace/extension curves or the retrace/retraction curves should be summarized.
        summary_func : function
            Function that will perform the summary. By default, this is a function that takes `z_piezo`, 
            `tm_defl`, and `sc_params` and returns the compliance matrix and R^2 matrix (how well each curve 
            was summarized). 
        **kwargs : dict
            Arguments that are passed to `summary_func`.

        Returns
        -------
        Default return values if `summary_func`=`get_comp_mat`.
        ndarray
            Compliance matrix. Shape should be square, with the size of the 
            sides being the square root of the number of force ramps.
            E.g. shape is (64, 64).
        ndarray
            R^2 matrix. See `comp` for shape.

        See Also
        --------
        pyMFD.summarize.get_comp_mat
        '''
        if which_dir == 'trace' or which_dir == 'extend' or which_dir == 0:
            which_dir = 0
        elif which_dir == 'retrace' or which_dir == 'retract' or which_dir == 1:
            which_dir = 1
        else:
            raise ValueError("which_dir must be one of: trace, extend, 0; or retrace, retract, 1")

        data = self.get_extend() if which_dir == 0 else self.get_retract()

        return summary_func(self.z_piezo, data, self.sc_params, **kwargs)

