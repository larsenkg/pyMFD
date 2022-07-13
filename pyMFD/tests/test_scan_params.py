from pyMFD.scan_params import get_scan_params
import os

def test_get_scan_params():
    '''Test loading scan parameters from JSON.'''
    sc_params = get_scan_params("data/examples/02041411.001.json")
    
    assert "name" in sc_params