"""Source/spec integration contract, not biological response calibration."""
import numpy as np
from sugarcode.modules.syndroid.core import mechanistic_cell_program


def test_seed_changes_stochastic_expression_output():
    low = mechanistic_cell_program(hours=1, seed=2)
    high = mechanistic_cell_program(hours=1, seed=3)
    assert low['stochastic_expression'] != high['stochastic_expression']
    assert low['stochastic_expression']['total_event_count'] > 0
    assert high['stochastic_expression']['total_event_count'] > 0


def test_stochastic_expression_is_coupled_to_cell_trajectory():
    low = mechanistic_cell_program(hours=1, seed=2)
    high = mechanistic_cell_program(hours=1, seed=3)
    # The function describes coupled expression/metabolism/allocation/feedback.
    # Its two genuinely different SSA outputs should have a path into the cell
    # trajectory, if this is an integrated stochastic cell model.
    assert any(not np.array_equal(low['states'][key], high['states'][key])
               for key in ('mRNA', 'protein', 'ATP', 'mass', 'circuit_output'))
