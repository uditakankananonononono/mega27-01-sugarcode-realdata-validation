from sugarcode.modules.vector_opt import engineer_capsid
from sugarcode.modules.cell_free_opt import kinetics
def test_receptor_identity_affects_model():
    a=engineer_capsid(target_receptor='AAVR',seed=1)
    b=engineer_capsid(target_receptor='heparan_sulfate',seed=1)
    assert a['variants']!=b['variants']
def test_initial_produced_mass_is_zero():
    r=kinetics({'yield_g_l':2},hours=0)
    assert r['trajectory'][0]['cumulative_g_l']==0 and r['final_g_l']==0
def test_seeded_positive_control():
    assert engineer_capsid(seed=1)==engineer_capsid(seed=1)
