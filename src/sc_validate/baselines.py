"""Naive reference baselines so every module result is compared against
something a skeptic would try first (mean predictor, simple physics proxy)."""
import csv, math, pathlib
from scipy.stats import pearsonr, spearmanr

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
# Kyte-Doolittle hydropathy
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
# residue volumes (A^3, Zamyatnin)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2,
           I=166.7, L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8,
           Y=193.6, V=140.0)


def load_skempi_subset(path=DATA / "skempi_subset_150.tsv"):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = r["mutation"]
            rows.append((m[0], m[-1], float(r["ddg_kcal"])))
    return rows


def skempi_baselines(rows):
    y = [d for _, _, d in rows]
    mean = sum(y) / len(y)
    rmse_mean = math.sqrt(sum((v - mean) ** 2 for v in y) / len(y))
    x_vol = [abs(VOL[w] - VOL[m]) for w, m, _ in rows]
    x_kd = [KD[w] - KD[m] for w, m, _ in rows]
    return {
        "n": len(y),
        "mean_ddg": round(mean, 4),
        "rmse_mean_predictor": round(rmse_mean, 4),
        "pearson_abs_volume_change": round(float(pearsonr(x_vol, y)[0]), 4),
        "spearman_abs_volume_change": round(float(spearmanr(x_vol, y)[0]), 4),
        "pearson_hydropathy_loss": round(float(pearsonr(x_kd, y)[0]), 4),
    }
