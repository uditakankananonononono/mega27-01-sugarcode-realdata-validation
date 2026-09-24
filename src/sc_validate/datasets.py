"""Loaders for held-out external datasets. All fixtures carry provenance."""
import csv, json, pathlib

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"


def clinvar_sample(path=None):
    """500-variant stratified ClinVar sample (250 pathogenic-class, 250 benign-class).

    Source: NCBI ClinVar variant_summary.txt.gz (FTP, retrieved 2026-09-24),
    GRCh38 SNVs with 'criteria provided' review status, reservoir-sampled.
    """
    p = pathlib.Path(path) if path else DATA / "clinvar_sample_500.tsv"
    rows = []
    for line in p.read_text().splitlines():
        f = line.split("\t")
        rows.append({
            "variation_id": f[0], "gene": f[1], "hgvs": f[2],
            "significance": f[3], "review": f[4], "chrom": f[5],
            "pos": f[6], "ref": f[7], "alt": f[8],
        })
    return rows


def skempi_single_point(path=None):
    """SKEMPI 2.0 single-point mutations with wt+mut affinities (n=4956)."""
    import math
    p = pathlib.Path(path) if path else DATA / "skempi_v2.csv"
    rows = list(csv.reader(open(p), delimiter=";"))
    h = rows[0]; idx = {k: i for i, k in enumerate(h)}
    out = []
    for r in rows[1:]:
        mut = r[idx["Mutation(s)_cleaned"]]
        if "," in mut or not mut:
            continue
        try:
            ka_wt = float(r[idx["Affinity_wt_parsed"]])
            ka_mut = float(r[idx["Affinity_mut_parsed"]])
        except ValueError:
            continue
        # ddG = R*T*ln(Kd_mut/Kd_wt) at 300K -> RT=0.596 kcal/mol
        out.append({
            "pdb": r[idx["#Pdb"]], "mutation": mut,
            "ddg_kcal": 0.596 * math.log(ka_mut / ka_wt),
        })
    return out


def bigg_e_coli_core(path=None):
    """BiGG e_coli_core model JSON (Orth et al. 2011 protocol model)."""
    p = pathlib.Path(path) if path else DATA / "e_coli_core.json"
    return json.loads(p.read_text())
