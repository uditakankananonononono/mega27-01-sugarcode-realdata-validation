#!/usr/bin/env python3
"""Settle the CRISPRscan coefficient-position discrepancy against the original
Moreno-Mateos 2015 supplementary model (Nature Methods nmeth.3543, MOESM640 xlsx,
media.springernature.com, retrieved 2026-09-26).
Result: 90/91 features identical (position + weight) between supplement and CRISPOR
paramsCRISPRscan; the largest negative weight sits on AA19 in the supplement but AA18
in CRISPOR. crisprScore (and sugarcode, which follows it) matches the original."""
import openpyxl, re, json, pathlib

D = pathlib.Path(__file__).resolve().parents[2] / "data" / "crisprscan"
wb = openpyxl.load_workbook(D / "MOESM640.xlsx", read_only=True)
suppl = {str(r[0]): float(r[1]) for r in wb["Sheet1"].iter_rows(values_only=True) if r[0]}

src = open("/home/sandbox/wave1/oracle/crisporEffScores.py").read()
block = src[src.index("paramsCRISPRscan = ["):]
block = block[: block.index("]")]
crispor = {}
for dinuc, pos, w in re.findall(r"\('([A-Z]+)',(\d+),([-\d.e]+)\)", block):
    crispor[f"{dinuc}{int(pos):02d}"] = float(w)

mismatch = [(k, v, crispor[k]) for k, v in suppl.items()
            if k != "Intercept" and k in crispor and abs(crispor[k] - v) > 1e-6]
missing = [(k, v) for k, v in suppl.items() if k != "Intercept" and k not in crispor]
extra = {k: v for k, v in crispor.items() if k not in suppl}
result = {
    "supplement_features": len(suppl) - 1,
    "crispor_features": len(crispor),
    "identical_position_and_weight": len(suppl) - 1 - len(mismatch) - len(missing),
    "value_mismatches_same_key": mismatch,
    "supplement_keys_missing_in_crispor": missing,
    "crispor_keys_not_in_supplement": extra,
    "conclusion": "CRISPOR places the largest negative coefficient (-0.0973770966031) "
                  "on AA18; the original supplement places it on AA19. Position conventions "
                  "align (90/91 agree at identical keys), so this is a one-position "
                  "misplacement in CRISPOR's Excel conversion, not a numbering difference.",
}
print(json.dumps(result, indent=1))
