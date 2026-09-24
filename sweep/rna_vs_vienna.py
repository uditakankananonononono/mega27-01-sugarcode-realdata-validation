"""rna_nussinov vs ViennaRNA 2.7 MFE on Rfam seed sequences with their
consensus structure projected as reference; base-pair F1 (Mathews 2004)."""
import json, os, sys, random, pathlib, re
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, os.environ.get("SUGARCODE_SRC", str(ROOT.parent / "sc-ai" / "src")))
import RNA
from sugarcode.modules.rna_nussinov import core as nu
OK = {("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
def read_sto(p):
    seqs, ss = {}, ""
    for l in open(p):
        if l.startswith("#=GC SS_cons"): ss += l.split()[-1]
        elif l.strip() and not l.startswith("#") and not l.startswith("//"):
            n, s = l.split()[:2]; seqs[n] = seqs.get(n, "") + s
    return seqs, ss
def cons_pairs(ss):
    st = {}; pairs = []
    br = {")": "(", "]": "[", "}": "{", ">": "<"}
    for i, c in enumerate(ss):
        if c in "([{<": st.setdefault(c, []).append(i)
        elif c in br and st.get(br[c]): pairs.append((st[br[c]].pop(), i))
    return pairs
def project(seq, ss):
    idx = {}; out = []
    for col, ch in enumerate(seq):
        if ch not in "-.~": idx[col] = len(out); out.append(ch.upper().replace("T", "U"))
    s = "".join(out)
    ref = {(idx[i], idx[j]) for i, j in cons_pairs(ss) if i in idx and j in idx and (s[idx[i]], s[idx[j]]) in OK}
    return s, ref
def f1(pred, ref):
    tp = len(pred & ref); p = tp / len(pred) if pred else 0; r = tp / len(ref) if ref else 0
    return 2 * p * r / (p + r) if p + r else 0.0
def db_pairs(db):
    st, out = [], set()
    for i, c in enumerate(db):
        if c == "(": st.append(i)
        elif c == ")": out.add((st.pop(), i))
    return out
res = {}
random.seed(7)
for fam in ("RF00005", "RF00001", "RF00010"):
    seqs, ss = read_sto(ROOT / "data" / "rfam" / f"{fam}.sto")
    items = [project(s, ss) for s in seqs.values()]
    items = [(s, r) for s, r in items if set(s) <= set("ACGU") and len(s) <= 400 and len(r) >= 5]
    items = random.sample(items, min(25, len(items)))
    fn, fv = [], []
    for s, ref in items:
        r = nu.fold(s, count_optimal=False)
        pn = set(tuple(p) for p in (r["pairs"] if isinstance(r, dict) else r.pairs))
        fn.append(f1(pn, ref)); fv.append(f1(db_pairs(RNA.fold(s)[0]), ref))
    res[fam] = {"n": len(items), "nussinov_f1": round(sum(fn) / len(fn), 4), "vienna_mfe_f1": round(sum(fv) / len(fv), 4),
                "mean_len": round(sum(len(s) for s, _ in items) / len(items), 1)}
    print(fam, res[fam], flush=True)
json.dump(res, open(ROOT / "benchmarks" / "sweep_rna.json", "w"), indent=1)
