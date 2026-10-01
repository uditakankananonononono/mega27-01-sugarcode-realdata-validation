"""Independent validation evidence expected to fail at product 1ccbd030."""
from sugarcode.modules.riboswitch import response_curve

def test_zero_basal_is_handled_without_unexpected_arithmetic_crash():
    # The public validator explicitly allows zero. A documented ValueError
    # would be acceptable if that contract changes; ZeroDivisionError is not.
    try:
        result = response_curve(1, basal=0, concentrations=[0, 1, 10])
    except ValueError:
        return
    assert result['curve'][1]['occupancy'] == .5
