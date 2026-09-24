"""mit_offtarget and crisprscan_score vs verbatim CRISPOR reference code on
held-out BRCA1 (NG_005905.2:60000-70000) guides."""
import csv, json, os, sys, re, random, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, os.environ.get("SUGARCODE_SRC", str(ROOT.parent / "sc-ai" / "src")))
sys.path.insert(0, str(ROOT / "sweep"))
from _crispor_mit_ref import calcHitScore
from sugarcode.modules.mit_offtarget import core as mit
from sugarcode.modules.crisprscan_score import core as cs
pairs = list(csv.DictReader(open(ROOT / "data" / "crispr_fresh_pairs.tsv"), delimiter="\t"))
d = []; exact = 0
for p in pairs:
    g, o = p["guide"][:20], p["offtarget"][:20]
    ref = calcHitScore(g, o) / 100.0
    mod = mit.score(g, o, "NGG").score
    d.append(abs(ref - mod)); exact += abs(ref - mod) < 1e-9
seq = "".join(l.strip() for l in open(ROOT / "data" / "brca1_10kb.fa") if not l.startswith(">")).upper()
ctx = [seq[m.start() - 6 - 20 + 0: m.start() - 20 + 20 + 3 + 6] for m in re.finditer(r"(?=[ACGT]GG)", seq)]
ctx = [c for c in ctx if len(c) == 35 and "N" not in c]
random.seed(7); ctx = random.sample(ctx, 300)
src = open(ROOT / "sweep" / "_crispor_effscores_ref.py").read()
ns = {}
i = src.index("paramsCRISPRscan"); j = src.index("def listToSvml")
exec(src[i:j], ns)
ref = ns["calcCrisprScanScores"](ctx)
mods = [cs.score(c) for c in ctx]
modv = [getattr(m, "score", None) for m in mods]
raw_ref = [r for r in ref]
agree = sum(int(r) == int(100 * v + 1e-9) for r, v in zip(ref, modv))
max_int_diff = max(abs(int(r) - int(100 * v + 1e-9)) for r, v in zip(ref, modv))
res = {"mit_pairs": len(pairs), "mit_exact": exact, "mit_max_abs_diff": max(d),
       "crisprscan_contexts": len(ctx), "crisprscan_int_agree": agree, "crisprscan_max_int_diff": max_int_diff,
       "crisprscan_example": [ctx[0], ref[0], modv[0]]}
json.dump(res, open(ROOT / "benchmarks" / "sweep_crispr.json", "w"), indent=1); print(res)
