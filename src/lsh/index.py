"""Step 5 of MMDS 3.4.3: find candidate pairs with the banding technique (3.4.1)."""

from collections import defaultdict
from collections.abc import Hashable
from itertools import combinations

import numpy as np


class LSHIndex[K: Hashable]:
    """Splits each signature into ``bands`` bands of ``rows`` rows and hashes each band.

    Two documents become a candidate pair when they land in the same bucket for
    at least one band. Bucketing by the band's exact values is the book's
    "hash to a large number of buckets", with no accidental collisions.
    """

    def __init__(self, bands: int, rows: int):
        self.bands = bands
        self.rows = rows
        self._buckets: defaultdict[tuple, list[K]] = defaultdict(list)

    def add(self, key: K, signature: np.ndarray) -> None:
        if signature.size != self.bands * self.rows:
            raise ValueError(
                f"signature has {signature.size} components, expected bands * rows = {self.bands * self.rows}"
            )
        for band, values in enumerate(signature.reshape(self.bands, self.rows)):
            self._buckets[band, tuple(values.tolist())].append(key)

    def candidate_pairs(self) -> set[tuple[K, K]]:
        return {pair for keys in self._buckets.values() for pair in combinations(keys, 2)}
