"""Held-out check: BiGG e_coli_core aerobic-glucose max growth must equal the
published 0.8739 h^-1 (Orth et al. 2011, PMID 21988831; BiGG benchmark)."""
from sc_validate.datasets import bigg_e_coli_core
from sc_validate.fba import solve_max_biomass


def test_e_coli_core_max_growth_matches_published():
    growth, flux, rxns = solve_max_biomass(bigg_e_coli_core())
    assert abs(growth - 0.873922) < 1e-4
    glc = flux[rxns.index("EX_glc__D_e")]
    assert abs(glc + 10.0) < 1e-6  # glucose-limited at 10 mmol/gDW/h
