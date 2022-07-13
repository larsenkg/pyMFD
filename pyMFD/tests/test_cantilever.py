import pytest
from pyMFD.FV import FV
from pyMFD.cantilever import get_cantilever_pos, get_cantilever_params, get_compliance_row, fit_compliance_linear, calc_modulus_offset, fit_compliance

def test_all():
    spm_file        = "data/examples/02041411.001"    # Example force-volume scan
    fv              = FV(spm_file)                    # Load force-volume scan
    (comp_mat, r2s) = fv.summarize()
    pos             = get_cantilever_pos(fv.get_pixel_size(), comp_mat.shape[0])

    (thick, width, igno, fixed, start, end, row, col_s, col_e) = get_cantilever_params(fv.sc_params, 0)

    comp_row           = get_compliance_row(comp_mat, row, rows_to_avg = 3)
    (slope, intercept) = fit_compliance_linear(pos[col_s:col_e], comp_row[col_s:col_e])
    (E_lin, off_lin)   = calc_modulus_offset(slope, intercept, width, thick)
    (E, offset, a)     = fit_compliance(pos[col_s:col_e], comp_row[col_s:col_e]**3, width, thick)

    assert pos.shape == (64,)
    assert thick == pytest.approx(1.6e-07)
    assert comp_row.shape == (64,)
    assert slope is not None
    assert E_lin is not None
    assert E is not None