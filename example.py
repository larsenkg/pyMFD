#from pyMFD.nanoscope import read_spm_data, read_spm_header, convert_spm_data, get_useful_params
from pyMFD.FV import FV
from pyMFD.summarize import comp_mat_inspector

import matplotlib.pyplot as plt

spm_file   = "data/examples/02041411.001"

fv = FV(spm_file)

print(fv.tm_defl.shape)

(comp_mat, r2s) = fv.summarize()
print(comp_mat.shape)

comp_mat_inspector(comp_mat, fv.z_piezo, fv.get_retract(), fv.sc_params, r2s_mat = r2s)
plt.show()

# fig, ax = plt.subplots()
# ax.pcolormesh(comp_mat, vmin=0, vmax=1)
# ax.invert_yaxis()
# plt.show()

# spm_params = read_spm_header(spm_file)
# params     = get_useful_params(spm_params)
# spm_data   = read_spm_data(spm_file, params)

# # Convert to metric units
# (z_piezo, tm_defl) = convert_spm_data(spm_data, params)



# print(params)