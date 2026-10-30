"""Step 4 of MMDS 3.4.3: choose bands b and rows r for a threshold t (3.4.2)."""


def candidate_probability(s: float, bands: int, rows: int) -> float:
    """Probability that two documents of Jaccard similarity ``s`` become a candidate pair."""
    return 1 - (1 - s**rows) ** bands


def threshold(bands: int, rows: int) -> float:
    """Similarity at which the S-curve is steepest: approximately (1/b)^(1/r)."""
    return (1 / bands) ** (1 / rows)


def choose_bands(n: int, t: float) -> tuple[int, int]:
    """Pick (bands, rows) with bands * rows == n whose threshold is closest to ``t``."""
    factorisations = [(n // r, r) for r in range(1, n + 1) if n % r == 0]
    return min(factorisations, key=lambda br: abs(threshold(*br) - t))
