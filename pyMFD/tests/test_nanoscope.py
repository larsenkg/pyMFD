import pytest
import struct
from pyMFD.nanoscope import get_fv_data, get_params, save_txt_data

### get_fv_data ###
# Tests get_fv_data, read_fv_data, convert_fv_data, get_params, read_fv_header, and convert_params

def test_get_fv_data():
    filename = "data/examples/02041411.001"
    params   = get_params(filename)
    data     = get_fv_data(filename, params)
    assert data is not None and params is not None

def test_get_fv_data_not_exist():
    with pytest.raises(FileNotFoundError) as e_info:
        filename = "data/examples/02041411.001"
        params   = get_params(filename)
        data     = get_fv_data("this_file_should_not_exist.fn238ty0asd9", params)

def test_get_fv_data_bad_filetype():
    with pytest.raises(struct.error) as e_info:
        filename = "data/examples/02041411.001"
        params   = get_params(filename)
        data     = get_fv_data("data/examples/02041411.001.json", params)

