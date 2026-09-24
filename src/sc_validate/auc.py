"""Rank-based AUC (Mann-Whitney) with no dependencies."""


def auc(scores_pos, scores_neg):
    """Probability a random positive outranks a random negative."""
    wins = 0.0
    for s in scores_pos:
        for t in scores_neg:
            wins += 1.0 if s > t else 0.5 if s == t else 0.0
    return wins / (len(scores_pos) * len(scores_neg))
