"""Binary-classification and agreement metrics with bootstrap CIs."""
import math, random


def confusion(y_true, y_pred):
    tp = sum(1 for t, p in zip(y_true, y_pred) if t and p)
    fp = sum(1 for t, p in zip(y_true, y_pred) if not t and p)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t and not p)
    tn = sum(1 for t, p in zip(y_true, y_pred) if not t and not p)
    return tp, fp, fn, tn


def binary_metrics(y_true, y_pred):
    tp, fp, fn, tn = confusion(y_true, y_pred)
    n = tp + fp + fn + tn
    return {
        "n": n, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "accuracy": (tp + tn) / n if n else float("nan"),
        "sensitivity": tp / (tp + fn) if tp + fn else float("nan"),
        "specificity": tn / (tn + fp) if tn + fp else float("nan"),
    }


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    return cov / math.sqrt(vx * vy)


def bootstrap_ci(stat, xs, ys=None, iters=1000, seed=42, alpha=0.05):
    rng = random.Random(seed)
    n = len(xs)
    vals = []
    for _ in range(iters):
        idx = [rng.randrange(n) for _ in range(n)]
        if ys is None:
            vals.append(stat([xs[i] for i in idx]))
        else:
            vals.append(stat([xs[i] for i in idx], [ys[i] for i in idx]))
    vals.sort()
    lo = vals[int(alpha / 2 * iters)]
    hi = vals[int((1 - alpha / 2) * iters)]
    return lo, hi
