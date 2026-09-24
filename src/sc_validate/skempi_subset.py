"""Select the SKEMPI subset used for the mutdock ddG check."""
import csv, pathlib


def pick_subset(path, n=150, seed=7):
    import random
    p = pathlib.Path(path)
    rows = list(csv.reader(open(p), delimiter=";"))
    h = rows[0]; idx = {k: i for i, k in enumerate(h)}
    cand = [r for r in rows[1:]
            if r[idx["Mutation(s)_cleaned"]] and "," not in r[idx["Mutation(s)_cleaned"]]
            and r[idx["Affinity_mut_parsed"]] and r[idx["Affinity_wt_parsed"]]]
    rng = random.Random(seed)
    rng.shuffle(cand)
    return cand[:n]
