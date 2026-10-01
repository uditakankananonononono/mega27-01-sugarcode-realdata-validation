"""Stress-conditioned survival contract, not measured apoptosis validation."""
import math
import numpy as np
from sugarcode.modules.bioprint_pro.core import oldroyd_b_extrusion


def test_pressure_increases_flow_and_reported_wall_stress():
    low = oldroyd_b_extrusion(20, .4, seed=5)
    high = oldroyd_b_extrusion(200, .4, seed=5)
    assert high['flow_rate_mm3_s'] > low['flow_rate_mm3_s']
    assert high['shear_stress_Pa'] > low['shear_stress_Pa']
    assert np.array_equal(low['cell_radial_position_mm'], high['cell_radial_position_mm'])


def test_pressure_stress_conditions_cell_viability():
    low = oldroyd_b_extrusion(20, .4, seed=5)
    high = oldroyd_b_extrusion(200, .4, seed=5)
    # Spec59 states V=exp(-lambda*sigma), sigma mechanical stress. Actual code
    # uses local shear RATE times passage time, cancelling extrusion velocity.
    # This assertion encodes that advertised stress-conditioned model contract,
    # not a measured real-cell dose response or requirement for every real ink.
    assert high['mean_viability'] < low['mean_viability']
    assert not math.isclose(high['mean_viability'], low['mean_viability'], abs_tol=1e-12)
