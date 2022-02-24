import numpy as np

from struct import unpack


FFL  = b'*Force file list'
CFIL = b'*Ciao force image list'
CIL  = b'*Ciao image list'
SL   = b'*Scanner list'
CSL  = b'*Ciao scan list'

def read_fv_header(filename: str) -> dict:
    '''
    Read the header information from a Bruker/Veeco Nanoscope v7.2 file. Returns a dictionary containing all of the lines
    from the header organized under the sections:
     - FFL  = b'*Force file list'
     - CFIL = b'*Ciao force image list'
     - CIL  = b'*Ciao image list'
     - SL   = b'*Scanner list'
     - CSL  = b'*Ciao scan list'

    Nanoscope header files are a mess. There will be different sections depending on the type of data in the file. For more
    information see Nanoscope User Guide and this informative forum post:
     - https://physics-astronomy-manuals.wwu.edu/Nanosocpe%207.3%20User%20Guide.pdf (broken link as of 2/23/2022)
     - http://nanoqam.ca/wiki/lib/exe/fetch.php?media=nanoscope_software_8.10_user_guide-d_004-1025-000_.pdf
     - http://nanoscaleworld.bruker-axs.com/nanoscaleworld/forums/p/538/1065.aspx

    In the file header some parameters start with '\@' instead of simply '\'. This is an indication to the software
    that the data that follows is intended for a CIAO parameter object. After the '@', you might see a number
    followed by a colon before the label. This number is what we call a “group number” and can generally be
    ignored.
    
    Further, after the label and its colon, you will see a single definition character of 'V', 'C', or 'S'.
     - V means _Value_ -- a parameter that contains a double and a unit of measure, and some scaling definitions.
     - C means _Scale_ -- a parameter that is simply a scaled version of another.
     - S means _Select_ -- a parameter that describes some selection that has been made

    '''

    # May be usefule:
    #  \*Ciao scan list\Scan Size: 15000 nm
    #  \*Ciao scan list\Samps/line: 64   # Force volume must be square.

    # Some code adapted from 'pySPM'. That library did not support force-volume
    # data, so I wrote this.
    # https://github.com/scholi/pySPM/blob/master/pySPM/Bruker.py
    params = {}
    key    = ""
    index  = 0
    count = 0
    
    with open(filename, 'rb') as file:
        while True:
            line = file.readline().rstrip().replace(b'\\', b"")
            if line[0] == 42: # Checking for asterisk, but line[0] returns integer
                if line == b'*File list end':
                    break
                if line not in params:
                    key         = line
                    params[key] = [{}]            
                    index       = 0
                else:
                    key   = line
                    index = len(params[key])
            else:
                args = line.split(b": ")
                if len(args) == 2:
                    params[key][index][args[0]] = args[1]
                else:
                    params[key][index][args[0]] = ""
            
            # TO DO: I should be able to remove this failsafe.
            # Don't read more than 1000 lines.
            count += count
            if count > 1000:
                break
            

    # Warn if version isn't 7.2
    sup_ver = b'0x07200000'
    if params[FFL][0][b'Version'] != sup_ver:
        print(f"Warning: Unsupported version detected. pyMFD only supports Nanoscope 7.2 ({sup_ver}), but version {params[FFL][0][b'Version']} detected.")
    # If CIL exists, then this is FV data. If it doesn't, this is a single force ramp curve.
    if CIL in params:
        params["is_single_curve"] = False
    else:
        params["is_single_curve"] =  True
    return params

def convert_params(old_params, custom_to_extract = []):
    '''
    CFIL
     - Data offset
     - Data length
     - Bytes/pixel
     - Samps/line
     - @4:Ramp size

    CSL
     - Samps/line
     - @2:TMDeflectionLimit

    SL
     - @Sens. Zsens
    '''

    from_value_f = lambda x: float(x.split()[-2])
    from_value_i = lambda x: int(x.split()[-2])

    to_extract = [
        # Section, Parameter Name         , New parameter name  , Function to convert from bytestring to desired type
            (CFIL, b"Data offset"         , "fv_data_offset"    , int),
            (CFIL, b"Data length"         , "fv_data_length"    , int),
            (CFIL, b"Bytes/pixel"         , "fv_bytes_per_pixel", int),
            (CFIL, b"Samps/line"          , "samples_per_ramp"  , from_value_i),
            (CFIL, b"@4:Ramp size"        , "ramp_size"         , from_value_f),
            (CSL,  b"Scan Size"           , "scan_size"         , lambda x: from_value_i(x)*1e-9), # convert from nm to m
            (CSL,  b"Samps/line"          , "ramps_per_line"    , int),
            (CSL,  b"@2:TMDeflectionLimit", "tm_volt_limit"     , from_value_f),
            (SL,   b"@Sens. Zsens"        , "piezo_nm_per_volt" , from_value_f),
            (None, "is_single_curve"      , "is_single_curve"   , lambda x: x)
    ]

    # Add any custom parameters to extract
    to_extract += custom_to_extract

    params = {}
    for (section, name, new_name, conv_func) in to_extract:
        if section is None:
            sec = old_params
        else:
            sec = old_params[section][0]

        params[new_name] = conv_func(sec[name])

    return params
    

def read_fv_data(filename: str, params: dict) -> np.ndarray:
    '''
    Read the force-volume or force-ramp data from a Nanoscope file. The data is converted from binary
    representation to a float64 representation of the the SPM data in ADC counts. Convert to volts
    using `convert_fv_data`.

    A force-volume scan contains three dimensions of data. For every point in a 2D array, two force-ramps are 
    recorded (one for extension towards the sample and one for retraction -- also called trace and retrace).

    The raw data should have a size equal to the number of points in the 2D array times the number of samples in 
    the force-ramp all times two (for extend and retract).

    For example, a 64x64 with 1024 samples per force-ramp will have a data length of:
     - 64^2 * 1024 * 2 = 8388608

     This length should be recorded in the header as `\*Ciao force image list\Data length` (keeping in mind the bytes/pixel).
    '''
    offset      = params["fv_data_offset"]
    data_length = params["fv_data_length"]
    bytes_px    = params["fv_bytes_per_pixel"]

    if params["is_single_curve"]:
        num_curves = 1
    else:
        num_curves = params["ramps_per_line"]**2
    
    # Two times the samples_per_ramp because their are trace and retrace (or extend and retract)
    ramp_length  = 2*params["samples_per_ramp"]
    
    fmt_chars  = {1: "b", 2: "h", 4: "i", 8: "q"}   # This is untested except for 2: "h"
    unpack_fmt = "<"                                # Little-endian
    unpack_fmt += str(ramp_length*num_curves)       # Number of bytes
    unpack_fmt += fmt_chars[bytes_px]               # Size of byte

    with open(filename, 'rb') as file:
        file.seek(offset)
        
        raw_data = file.read(data_length)

    return np.array(unpack(unpack_fmt, raw_data), dtype='float64')
    
    
def convert_fv_data(data: np.ndarray, params: str) -> tuple:
    '''
    Convert from ADC counts to volts. Returns the piezo ramp deflection `z_piezo` and the 
    force-volume TM deflection data in volts in a tuple: (z_piezo, tm_defl).
    '''
    z_sens     = params["piezo_nm_per_volt"]
    ramp_size  = params["ramp_size"]
    tm_limit   = params["tm_volt_limit"]
    ramp_len   = params["samples_per_ramp"]
    bits       = params["fv_bytes_per_pixel"]*8   # Should be 16
    
    if params["is_single_curve"]:
        num_curves = 1
    else:
        num_curves = params["ramps_per_line"]**2
            
    z_size     = z_sens * ramp_size / ramp_len * (ramp_len - 1)              # nanometers
    tm_mod     = tm_limit/2**bits                                            # volts / bit
    z_piezo    = np.linspace(0, z_size, ramp_len)                            # nanometers
    tm_defl    = data.reshape((ramp_len, 2, num_curves), order='F')          # bits (ADC counts)
    tm_defl   *= tm_mod                                                      # volts
 
    # Returns the ramp piezo deflection `z_piezo` in nanometers and `tm_defl` in volts.
    return (z_piezo, tm_defl)

def get_fv_data(filename: str, params: dict) -> tuple:
    '''Get the `z_piezo` deflection ramp. `params` should be the converted, generalized parameter dictionary.'''
    data  = read_fv_data(filename, params)

    # Convert to metric units
    return convert_fv_data(data, params)

def get_params(filename: str) -> dict:
    all_fv_params = read_fv_header(filename )
    fv_params     = convert_params(all_fv_params)

    return fv_params

def save_txt_data(data, filename):
    '''
    Save the converted data to an ASCII file using the same format as exports from Nanoscope Analysis 2.0.
    '''
    header = "Calc_Ramp_Ex_nm\tCalc_Ramp_Rt_nm\tDefl_mV_Ex\tDefl_mV_Rt\tpN Not Available\tpN Not Available\t"
    np.savetxt(filename, data, delimiter='\t', fmt='%1.6e', header=header, comments='')
