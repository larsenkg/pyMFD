#from pyMFD.nanoscope import read_spm_data, read_spm_header, convert_spm_data, get_useful_params
from pyMFD.FV import FV
from pyMFD.summarize import comp_mat_inspector

import matplotlib.pyplot as plt

spm_file = "data/examples/02041411.001"    # Example force-volume scan
fv       = FV(spm_file)                    # Load force-volume scan

print(fv.tm_defl.shape)

(comp_mat, r2s) = fv.summarize()
print(comp_mat.shape)

# Interactive compliance map inspector
# Use to mouse to select pixels in the (left) compliance map.
# The raw force-deflection data is shown in the center plot.
# The R^2 map (how well the force-deflection data was fit) is shown in the right map.
comp_mat_inspector(comp_mat, fv.z_piezo, fv.get_retract(), fv.sc_params)#, r2s_mat = r2s)
plt.show()
