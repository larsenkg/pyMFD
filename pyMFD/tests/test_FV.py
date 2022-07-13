import pytest
from pyMFD.FV import FV


def test_FV():
    spm_file  = "data/examples/02041411.001" 
    fv        = FV(spm_file)
    (comp, _) = fv.summarize()
    px_size   = fv.get_pixel_size()

    assert fv.tm_defl.shape == (1024, 2, 4096)
    assert comp.shape       == (64, 64)
    assert px_size          == pytest.approx(2.34375e-07)