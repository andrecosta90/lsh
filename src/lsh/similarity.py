"""Similarity measures used in steps 6 and 7 of MMDS 3.4.3."""

from collections.abc import Set

import numpy as np


def jaccard(a: Set, b: Set) -> float:
    """|A ∩ B| / |A ∪ B| (MMDS 3.1.1)."""
    return len(a & b) / len(a | b)


def signature_similarity(x: np.ndarray, y: np.ndarray) -> float:
    """Fraction of components in which two minhash signatures agree (MMDS 3.3.3)."""
    return float(np.mean(x == y))
