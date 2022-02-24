from pyMFD.nanoscope   import get_fv_data, get_params
from pyMFD.scan_params import get_scan_params

class FV:
    '''
    This class represents a single force-volume scan. It contains the relevant scan parameters and force-volume data.
    '''

    def __init__(self, fv_filename):
        self.fv_filename             = fv_filename
        self.fv_params               = get_params(self.fv_filename)
        (self.z_piezo, self.tm_defl) = get_fv_data(self.fv_filename, self.fv_params)
        self.sp_params               = get_scan_params(self.fv_filename + ".json")        
        self.pixel_size              = self.get_pixel_size()

    def get_pixel_size(self, scan_size=None, scan_points=None):
        if scan_size is None:
            scan_size   = self.fv_params["scan_size"]

        if scan_points is None:
            scan_points = self.fv_params["ramps_per_line"]

        return scan_size / scan_points

    def summarize(summary_func):
        '''Create a 2D representation of the force-volume data.'''

    