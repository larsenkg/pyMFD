#from pyMFD.nanoscope import read_spm_data, read_spm_header, convert_spm_data, get_useful_params
from pyMFD.FV import FV
from pyMFD.summarize import comp_mat_inspector
from pyMFD.cantilever import (
    get_cantilever_pos, 
    get_cantilever_params, 
    get_compliance_row, 
    fit_compliance_linear, 
    calc_modulus_offset, 
    standardize_and_fit
)
import matplotlib.pyplot as plt

use_inspector = True

spm_file = "data/examples/02041411.001"    # Example force-volume scan
fv       = FV(spm_file)                    # Load force-volume scan

print(fv.tm_defl.shape)

(comp_mat, r2s) = fv.summarize()
print(comp_mat.shape)

# Interactive compliance map inspector
# Use mouse to select pixels in the (left) compliance map.
# The raw force-deflection data is shown in the center plot.
# The R^2 map (how well the force-deflection data was fit) is shown in the right map.
if use_inspector:
    comp_mat_inspector(comp_mat, fv.z_piezo, fv.get_retract(), fv.sc_params, r2s_mat = r2s)
    plt.show()

###########################
# Find cantilever modulus #
###########################
cant_num    = 0 # Cantilever number
rows_to_avg = 3 # Number of rows to average around center line of cantilever
pos         = get_cantilever_pos(fv.get_pixel_size(), comp_mat.shape[0])

(thick, width, igno, fixed, start, end, row, col_s, col_e) = get_cantilever_params(fv.sc_params, cant_num)

comp_row           = get_compliance_row(comp_mat, row, rows_to_avg = rows_to_avg)
(slope, intercept) = fit_compliance_linear(pos[col_s:col_e], comp_row[col_s:col_e])
(E_lin, off_lin)   = calc_modulus_offset(slope, intercept, width, thick)
(E, offset, a)     = standardize_and_fit(pos[col_s:col_e], comp_row[col_s:col_e]**3, width, thick)

print(f"Width: {width*1e9:.2f} nm")
print(f"Thick: {thick*1e9:.2f} nm")

print("---- Cubic fit ----")
print(f"a: {a:.2g}")
print(f"Young's modulus: {E/1e9:.2f} GPa")
print(f"Offset: {offset*1e6:.2f} µm")

print("---- Linearized fit ----")
print(f"Slope: {slope:.2f}")
print(f"Y-int: {intercept:.2f}")
print(f"Young's modulus: {E_lin/1e9:.2f} GPa")
print(f"Offset: {off_lin*1e6:.2f} µm")
