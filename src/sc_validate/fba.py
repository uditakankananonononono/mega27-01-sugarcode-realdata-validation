"""Independent FBA reference (HiGHS via scipy) used to check virtual_cell parity."""
import numpy as np
from scipy.optimize import linprog


def solve_max_biomass(model, objective="BIOMASS_Ecoli_core_w_GAM"):
    mets = [x["id"] for x in model["metabolites"]]
    rxns = [x["id"] for x in model["reactions"]]
    mi = {k: i for i, k in enumerate(mets)}
    S = np.zeros((len(mets), len(rxns)))
    lb, ub = [], []
    for j, r in enumerate(model["reactions"]):
        for k, v in r["metabolites"].items():
            S[mi[k], j] = v
        lb.append(r["lower_bound"])
        ub.append(r["upper_bound"])
    c = np.zeros(len(rxns))
    c[rxns.index(objective)] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(mets)),
                  bounds=list(zip(lb, ub)), method="highs")
    assert res.success, res.message
    return -res.fun, res.x, rxns
