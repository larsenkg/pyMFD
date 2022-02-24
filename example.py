#from pyMFD.nanoscope import read_spm_data, read_spm_header, convert_spm_data, get_useful_params
from pyMFD.FV import FV

spm_file   = "data/examples/02041411.001"

fv = FV(spm_file)

print(fv.tm_defl.shape)
# spm_params = read_spm_header(spm_file)
# params     = get_useful_params(spm_params)
# spm_data   = read_spm_data(spm_file, params)

# # Convert to metric units
# (z_piezo, tm_defl) = convert_spm_data(spm_data, params)



# print(params)