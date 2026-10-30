"""Step 3 of MMDS 3.4.3: compress shingle sets into minhash signatures (3.3)."""

import zlib
from collections.abc import Iterable

import numpy as np

# Mersenne prime 2^31 - 1: the "c" in h(x) = (a*x + b) mod c (MMDS 3.3.5).
# Small enough that a*x + b never overflows uint64.
_PRIME = (1 << 31) - 1


class MinHasher:
    """Computes length-``n`` minhash signatures with ``n`` random hash functions.

    MMDS 3.3.5 builds signatures by scanning the rows of the characteristic
    matrix (one row per shingle, sorted) and keeping, for each hash function,
    the minimum hash of the rows where the column has a 1. That minimum depends
    only on the set of shingles in the column, so we compute it directly from
    the set and never materialise the (huge, sparse) matrix.
    """

    def __init__(self, n: int, seed: int | None = None):
        rng = np.random.default_rng(seed)
        self._a = rng.integers(1, _PRIME, size=(n, 1), dtype=np.uint64)
        self._b = rng.integers(0, _PRIME, size=(n, 1), dtype=np.uint64)

    def signature(self, shingles: Iterable[str]) -> np.ndarray:
        rows = np.fromiter((_row_number(s) for s in shingles), dtype=np.uint64)
        if rows.size == 0:
            raise ValueError("cannot compute the signature of an empty set")
        return ((self._a * rows + self._b) % _PRIME).min(axis=1)


def _row_number(shingle: str) -> int:
    """Hash a shingle to 4 bytes (MMDS 3.2.3), used as its row in the matrix."""
    return zlib.crc32(shingle.encode()) % _PRIME
