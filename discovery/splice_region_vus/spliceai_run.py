"""SpliceAI (Jaganathan 2019, official pip weights) delta scores for the
same CV variants, strand-oriented, D=50 (max over |pos|<=50 of the four
delta types), on a seeded stratified subset. Incremental output."""
import csv, sys, os, random, numpy as np
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import twobitreader
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import tensorflow as tf
tf.config.threading.set_inter_op_parallelism_threads(2); tf.config.threading.set_intra_op_parallelism_threads(2)
models = [load_model(resource_filename('spliceai', f'models/spliceai{i}.h5'), compile=False) for i in range(1, 6)]
RC = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(RC)[::-1]
g = twobitreader.TwoBitFile("hg38.2bit")
D = 50; W = 5000 + D
loc = {}
for r in csv.reader(open("splice_region_snvs.tsv"), delimiter="\t"):
    loc[r[0]] = (r[5], r[6], r[7], r[8])
win = {r['vid']: r for r in csv.DictReader(open('windows.tsv'), delimiter='\t')}
ids = [l.strip() for l in open(sys.argv[1])]
done = set()
if os.path.exists("spliceai_scores.tsv"):
    done = {l.split("\t")[0] for l in open("spliceai_scores.tsv")}
out = open("spliceai_scores.tsv", "a")
for n, vid in enumerate(ids):
    if vid in done: continue
    chrom, pos, ref, alt = loc[vid]; ch = "chr" + chrom; p = int(pos) - 1
    seq = g[ch][p - W:p + W + 1].upper()
    if len(seq) != 2 * W + 1 or seq[W] != ref: continue
    alts = seq[:W] + alt + seq[W + 1:]
    if win[vid]['strand'] == '-': seq, alts = rc(seq), rc(alts)
    x = np.stack([one_hot_encode(seq), one_hot_encode(alts)])
    y = np.mean([m.predict(x, verbose=0) for m in models], axis=0)  # (2, 2D+1, 3)
    r_, a_ = y[0], y[1]
    ag = (a_[:, 1] - r_[:, 1]).max(); al = (r_[:, 1] - a_[:, 1]).max()
    dg = (a_[:, 2] - r_[:, 2]).max(); dl = (r_[:, 2] - a_[:, 2]).max()
    out.write(f"{vid}\t{ag:.4f}\t{al:.4f}\t{dg:.4f}\t{dl:.4f}\t{max(ag,al,dg,dl):.4f}\n"); out.flush()
    if n % 50 == 0: print(n, flush=True)
