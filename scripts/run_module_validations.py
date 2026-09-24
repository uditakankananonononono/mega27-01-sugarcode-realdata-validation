#!/usr/bin/env python3
"""Run sugarcode module validations against frozen held-out fixtures.

Checks (frozen spec, 2026-09-24):
 1. crispr_opt CFD vs CRISPOR-extracted oracle, 600 fresh pairs
 2. crispr_opt RS1 vs 150 frozen CRISPOR percentiles
 3. openclinvar (+rarenet note) vs 500 stratified ClinVar labels
 4. deepsplice PWM AUC on held-out RefSeqGene junctions
 5. mutdock ddG vs SKEMPI subset (150) + naive baselines
 6. virtual_cell FBA parity vs 0.873922 (e_coli_core max biomass)
 7. profile_hmm globin-family AUC (Pfam PF00042 HMMER3 reference)
 8. codon CAI vs Kazusa published table
 9. acmg_bayesian vs Tavtigian worked examples

Writes benchmarks/results.json. sugarcode source root comes from
SUGARCODE_SRC env (e.g. /path/to/sugarcode-ai/src).
"""
import json, math, os, pathlib, random, re, sys, time, urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
sys.path.insert(0, os.environ.get("SUGARCODE_SRC", "/tmp/sc-ai/src"))
sys.path.insert(0, str(ROOT / "src"))

results = {"generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()), "checks": {}}

def record(name, status, metrics=None, note=""):
    results["checks"][name] = {"status": status, "metrics": metrics or {}, "note": note}
    print(f"[{status:>14}] {name}: {note[:120]}")

# ---------------------------------------------------------------- 1. CFD 600
from sugarcode.modules.crispr_opt import core as cr
from sc_validate import cfd_oracle

rows = [l.split("\t") for l in (DATA / "crispr_fresh_pairs.tsv").read_text().splitlines()[1:] if l.strip()]
agree = err = 0; maxdiff = 0.0
for g, o, _nmm in rows:
    try:
        m = cr.cfd_score(g[:20], o[:20], pam2=o[-2:])
        ref = cfd_oracle.cfd_score(g, o)
    except Exception:
        err += 1; continue
    d = abs(m - ref); maxdiff = max(maxdiff, d)
    if d <= 5e-7: agree += 1
record("cfd600", "pass" if agree == len(rows) else "fail",
       {"pairs": len(rows), "exact_agree": agree, "errors": err, "max_abs_diff": maxdiff},
       f"{agree}/{len(rows)} exact (tol 5e-7), max |diff| {maxdiff:.2e}")

# ---------------------------------------------------------------- 2. RS1 150
rs = [l.split("\t") for l in (DATA / "doench_rs1_oracle.tsv").read_text().splitlines()[1:] if l.strip()]
mods = [cr.doench2014_ontarget(ctx) for ctx, _ in rs]
refs = [float(p) for _, p in rs]

def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i]); r = [0.0] * len(xs); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]: j += 1
        for k in range(i, j + 1): r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r

def pearson(a, b):
    n = len(a); ma, mb = sum(a)/n, sum(b)/n
    cov = sum((x-ma)*(y-mb) for x, y in zip(a, b))
    return cov / math.sqrt(sum((x-ma)**2 for x in a) * sum((y-mb)**2 for y in b))

sp = pearson(ranks(mods), ranks(refs))
record("rs1_150", "pass" if sp > 0.5 else "fail",
       {"contexts": len(rs), "spearman": round(sp, 4),
        "module_range": [round(min(mods), 2), round(max(mods), 2)]},
       f"Spearman(module raw RS1, frozen CRISPOR percentile) = {sp:.3f}")

# ---------------------------------------------------------- 3. ClinVar 500
from sugarcode.modules.openclinvar import core as oc
from sc_validate import clinvar_runner

def classify(r):
    try:
        res = oc.interpret_variant(r["gene"], r["hgvs"])
    except Exception:
        return None
    cls = str(res.get("classification", "")).lower().replace("_", " ")
    if not cls or "uncertain" in cls or "vus" in cls: return None
    if "pathogenic" in cls: return "P"
    if "benign" in cls: return "B"
    return None

try:
    m = clinvar_runner.evaluate(classify)
except ZeroDivisionError:
    m = {"n": 0, "abstained": 500, "accuracy": None}
# Diagnosis (manual 8-variant live pilot, 2026-09-24): the live ClinVar query
# format fails to resolve the sample's HGVS strings, so enrichment never engages.
record("clinvar500", "fail" if (m.get("accuracy") or 0) < 0.5 else "pass", m,
       (f"acc {m['accuracy']:.3f} sens {m['sensitivity']:.3f} spec {m['specificity']:.3f}, " if m.get("accuracy") is not None else "accuracy undefined: ")
       + f"abstained {m['abstained']}/500 offline. Diagnosis: the consequence parser misses "
       "protein-style HGVS ('p.Arg518Ser' -> consequence 'unknown', neutral weight 0.35) and an "
       "8-variant live pilot (4 pathogenic, 4 benign) shows the ClinVar query format fails to resolve the sample's HGVS strings, so "
       "enrichment never engages; the module abstains on everything. Safety-conservative but "
       "degenerate as a classifier. rarenet is symptom-driven and emits no variant label: not scored.")

# -------------------------------------------- 4. deepsplice junction AUC
from sugarcode.modules.deepsplice import core as ds
from sc_validate import harvest_junctions as hj
from sc_validate import auc as svauc

TRAINING = set("""NG_012772 NG_017006 NG_007466 NG_008481 NG_016323 NG_033079
NG_009830 NG_012386 NG_034227 NG_054724 NG_005895 NG_008617 NG_005905 NG_007109
NG_017013 NG_007110 NG_007111 NG_008466 NG_016465 NG_008690 NG_009009 NG_009018
NG_009060 NG_059281 NG_007107 NG_007406 NG_008150 NG_011906 NG_008189 NG_016178
NG_007884 NG_008358 NG_011403 NG_012232 NG_008212 NG_011673""".split())
CANDIDATES = ["NG_008021.1", "NG_007489.1", "NG_029934.1", "NG_007485.1",
              "NG_009974.1", "NG_013013.2", "NG_008245.1", "NG_016442.1"]

def efetch_gb(acc):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
           f"{acc}&rettype=gb&retmode=text")
    return urllib.request.urlopen(url, timeout=30).read().decode()

def parse_genbank(text):
    m = re.search(r"^ORIGIN\s*\n(.*?)^//", text, re.S | re.M)
    seq = re.sub(r"[^acgtACGT]", "", m.group(1)).upper()
    cds = re.search(r"^     CDS\s+(join\(.*?\))", text, re.S | re.M)
    if not cds or "complement" in cds.group(1):
        return seq, []
    return seq, hj.parse_cds_exons("CDS " + cds.group(1))

try:
    don_pos, acc_pos, used = [], [], []
    for acc in CANDIDATES:
        if acc.split(".")[0] in TRAINING: continue
        try: text = efetch_gb(acc)
        except Exception: continue
        if "RefSeqGene" not in text: continue
        seq, exons = parse_genbank(text)
        js = hj.junction_windows(seq, exons)
        if len(js) < 6: continue
        for j in js:
            if j.kind == "donor" and len(j.window) == 9: don_pos.append(j.window)
            elif j.kind == "acceptor" and len(j.window) == 23: acc_pos.append(j.window[5:20])
        used.append(acc)
        if len(used) >= 3: break
    if not used: raise RuntimeError("no held-out RefSeqGene record fetched")
    rng = random.Random(7)
    def shuf(ws):
        out = []
        for w in ws:
            mid = list(w); rng.shuffle(mid); out.append("".join(mid))
        return out
    auc_d = svauc.auc([ds.score_donor(w) for w in don_pos], [ds.score_donor(w) for w in shuf(don_pos)])
    auc_a = svauc.auc([ds.score_acceptor(w) for w in acc_pos], [ds.score_acceptor(w) for w in shuf(acc_pos)])
    record("deepsplice_junction_auc", "pass" if auc_d > 0.8 and auc_a > 0.8 else "fail",
           {"accessions": used, "donor_junctions": len(don_pos), "acceptor_junctions": len(acc_pos),
            "donor_auc": round(auc_d, 4), "acceptor_auc": round(auc_a, 4)},
           f"held-out RefSeqGene {used} (disjoint from module PROVENANCE list); donor AUC {auc_d:.3f}, "
           f"acceptor AUC {auc_a:.3f} vs shuffled-window negatives (seed 7); plus-strand CDS only")
except Exception as e:
    record("deepsplice_junction_auc", "blocked", {}, f"fetch/score failed: {e}")

# ------------------------------------- 5. mutdock SKEMPI ddG (+ baselines)
from sc_validate.baselines import load_skempi_subset, skempi_baselines
bl = skempi_baselines(load_skempi_subset())
record("mutdock_skempi_ddg", "not_applicable", {"skempi_baselines": bl},
       "mutdock's ddG is a ligand-pocket class heuristic: mutation_effect(pocket_residues, smiles, "
       "position, mutant) mixes a pocket-ligand binding delta with a residue-class penalty. SKEMPI "
       "rows are protein-protein mutations with neither a ligand SMILES nor a pocket string, so no "
       "honest adapter exists; scoring would be theatre. Documented gap, contextualized by the naive "
       "baselines (mean-predictor RMSE, volume-change and hydropathy correlations) in metrics.")

# ------------------------------------------------------- 6. virtual_cell FBA
from sugarcode.modules.virtual_cell import core as vc
from sc_validate import fba as svfba
from sc_validate import datasets as svd
import numpy as np

mj = svd.bigg_e_coli_core()
mets = [x["id"] for x in mj["metabolites"]]; rxns = [x["id"] for x in mj["reactions"]]
mi = {k: i for i, k in enumerate(mets)}
S = np.zeros((len(mets), len(rxns))); lb = []; ub = []
for j, r in enumerate(mj["reactions"]):
    for k, v in r["metabolites"].items(): S[mi[k], j] = v
    lb.append(r["lower_bound"]); ub.append(r["upper_bound"])
mm = vc.MetabolicModel(mets, rxns, S, np.array(lb), np.array(ub),
                       objective="BIOMASS_Ecoli_core_w_GAM")
out = vc.fba(mm)
ref_obj, _, _ = svfba.solve_max_biomass(mj)
HEADLINE = 0.873922
ok = (out["status"] == "optimal" and abs(out["objective"] - ref_obj) < 1e-6
      and abs(ref_obj - HEADLINE) < 1e-6)
record("virtual_cell_fba_parity", "pass" if ok else "fail",
       {"module_objective": out["objective"], "independent_objective": round(ref_obj, 6),
        "headline": HEADLINE, "status": out["status"]},
       f"module {out['objective']} vs independent HiGHS {ref_obj:.6f} vs published headline {HEADLINE}")

# -------------------------------------- 7. globin AUC (Pfam/HMMER reference)
from sc_validate.pfam_reference import pfam_globin_reference
pf = pfam_globin_reference()
record("profile_hmm_globin_auc", "reference_only", pf,
       f"Pfam PF00042.29 (HMMER3) separates the fixtures perfectly (AUC {pf['auc']}), so the task "
       "is solvable and the fixture labels are sound. The module's profile_hmm still cannot be "
       "scored: it trains from an MSA and no globin alignment ships; AUC(module) deferred until "
       "the Clustal Omega alignment lands. Reference established, module verdict pending.")

# ------------------------------------------------------------- 8. CAI/Kazusa
from sugarcode.bio import codon as codonlib
pub = codonlib.load_published_table("e_coli_316407")
GC = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
      "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
      "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
      "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
      "TAT":"Y","TAC":"Y","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N",
      "AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C",
      "TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R",
      "AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
def indep_cai(seq):
    ws = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]; aa = GC.get(c)
        if aa in (None, "M", "W"): continue
        ws.append(pub[c] / max(pub[s] for s in GC if GC[s] == aa))
    return math.exp(sum(math.log(w) for w in ws) / len(ws))
rng = random.Random(11)
codons = [c for c, a in GC.items() if a and a not in ("M", "W")]
diffs = []
for _ in range(40):
    seq = "".join(rng.choice(codons) for _ in range(rng.randrange(60, 300)))
    diffs.append(abs(codonlib.cai(seq, pub) - indep_cai(seq)))
mx = max(diffs)
record("codon_cai_kazusa", "pass" if mx < 1e-6 else "fail", {"sequences": 40, "max_abs_diff": mx},
       f"module codon.cai vs independent Kazusa geometric-mean recompute, 40 random sequences, "
       f"max |diff| {mx:.2e}")

# ------------------------------------------------- 9. acmg_bayesian Tavtigian
from sugarcode.modules.acmg_bayesian import core as ac
OPVS = 350.0; PRIOR = 0.10
def indep_posterior(points, prior=PRIOR):
    odds = OPVS ** (points / 8) * (prior / (1 - prior))
    return odds / (1 + odds)
pmax = max(abs(ac.posterior_from_points(p) - indep_posterior(p)) for p in range(-16, 17))
bands = {"PATHOGENIC": [10, 12, 16], "LIKELY_PATHOGENIC": [6, 9], "VUS": [0, 5],
         "LIKELY_BENIGN": [-1, -6], "BENIGN": [-7, -12]}
band_ok = all(ac.classify_points(p).lower().replace(" ", "_") in (label.lower(), "uncertain_significance" if label == "VUS" else label.lower())
              for label, pts in bands.items() for p in pts)
worked = ac.posterior_from_points(6)
ok = pmax < 1e-12 and band_ok
record("acmg_tavtigian", "pass" if ok else "fail",
       {"max_posterior_diff": pmax, "bands_match_readme": band_ok,
        "six_point_posterior": round(worked, 4)},
       f"posterior matches OP=350^(pts/8) to {pmax:.1e} over -16..16; band labels "
       f"{'match' if band_ok else 'MISMATCH'}; 6 pts at prior 0.10 -> {worked:.3f} "
       "(Tavtigian 2018 prints this identity as 0.900)")

# ------------------------------------------------------------------ write
out = ROOT / "benchmarks"; out.mkdir(exist_ok=True)
(out / "results.json").write_text(json.dumps(results, indent=2))
n = {s: sum(1 for c in results["checks"].values() if c["status"] == s)
     for s in ("pass", "fail", "blocked", "not_applicable", "reference_only")}
print(f"\nwrote benchmarks/results.json: {n}")
