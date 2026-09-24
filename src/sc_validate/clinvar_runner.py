"""Run a variant classifier over the held-out ClinVar sample and score it."""
import csv, pathlib
from .metrics import binary_metrics, bootstrap_ci

PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}


def evaluate(classify_fn, sample_path=None, limit=None):
    """classify_fn(dict) -> 'P' | 'B' | None (None = module abstains)."""
    p = sample_path or (pathlib.Path(__file__).resolve().parents[2]
                        / "data" / "clinvar_sample_500.tsv")
    y_true, y_pred, abstained = [], [], 0
    rows = list(csv.DictReader(open(p), delimiter="\t",
                               fieldnames=["variation_id", "gene", "hgvs",
                                           "significance", "review", "chrom",
                                           "pos", "ref", "alt"]))
    if limit:
        rows = rows[:limit]
    for r in rows:
        pred = classify_fn(r)
        if pred is None:
            abstained += 1
            continue
        y_true.append(r["significance"] in PATHOGENIC)
        y_pred.append(pred == "P")
    m = binary_metrics(y_true, y_pred)
    lo, hi = bootstrap_ci(lambda a: sum(a) / len(a),
                          [t == p for t, p in zip(y_true, y_pred)])
    m["abstained"] = abstained
    m["accuracy_ci95"] = (round(lo, 4), round(hi, 4))
    return m
