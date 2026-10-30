"""Step 4 of MMDS 3.4.3: choose bands b and rows r for a threshold t (3.4.2)."""


def candidate_probability(s: float, bands: int, rows: int) -> float:
    """Probability that two documents of Jaccard similarity ``s`` become a candidate pair."""
    return 1 - (1 - s**rows) ** bands


def threshold(bands: int, rows: int) -> float:
    """Similarity at which the S-curve is steepest: approximately (1/b)^(1/r)."""
    return (1 / bands) ** (1 / rows)


def choose_bands(n: int, t: float, *, min_recall: float | None = None) -> tuple[int, int]:
    """Pick (bands, rows) with bands * rows == n.

    By default, the factorisation whose threshold is closest to ``t``. To avoid
    false negatives, MMDS suggests a threshold lower than ``t``: ``min_recall``
    says how much lower, by picking the highest threshold for which documents
    of similarity ``t`` still become candidates with probability ``min_recall``.
    """
    factorisations = [(n // r, r) for r in range(1, n + 1) if n % r == 0]
    if min_recall is None:
        return min(factorisations, key=lambda br: abs(threshold(*br) - t))

    reaching = [br for br in factorisations if candidate_probability(t, *br) >= min_recall]
    if not reaching:
        raise ValueError(f"no bands * rows = {n} reaches min_recall={min_recall} at t={t}; increase n")
    return max(reaching, key=lambda br: threshold(*br))
