"""MMDS 3.4.3, "Combining the Techniques": find similar documents end to end."""

from collections.abc import Hashable, Mapping
from dataclasses import dataclass

from lsh.banding import choose_bands
from lsh.index import LSHIndex
from lsh.minhash import MinHasher
from lsh.shingling import shingles
from lsh.similarity import jaccard, signature_similarity


@dataclass(frozen=True)
class SimilarPair[K: Hashable]:
    a: K
    b: K
    similarity: float


def find_similar[K: Hashable](
    docs: Mapping[K, str],
    *,
    k: int,
    n: int,
    t: float,
    seed: int | None = None,
    verify: bool = True,
    lowercase: bool = False,
    min_recall: float | None = None,
) -> list[SimilarPair[K]]:
    """Return the pairs of documents with similarity of at least ``t``, most similar first.

    With ``verify`` (step 7) similarity is the exact Jaccard similarity of the
    shingle sets; otherwise it is the minhash estimate from step 6.
    ``min_recall`` trades false positives for fewer false negatives when
    choosing bands (see ``choose_bands``).
    """
    shingle_sets = {key: shingles(text, k, lowercase=lowercase) for key, text in docs.items()}  # step 1

    minhasher = MinHasher(n, seed)  # steps 2-3
    signatures = {key: minhasher.signature(s) for key, s in shingle_sets.items()}

    index = LSHIndex(*choose_bands(n, t, min_recall=min_recall))  # step 4
    for key, signature in signatures.items():
        index.add(key, signature)

    def similarity(a: K, b: K) -> float:
        if verify:  # step 7
            return jaccard(shingle_sets[a], shingle_sets[b])
        return signature_similarity(signatures[a], signatures[b])

    pairs = [
        SimilarPair(a, b, similarity(a, b))
        for a, b in index.candidate_pairs()  # step 5
        if signature_similarity(signatures[a], signatures[b]) >= t  # step 6
    ]
    return sorted((p for p in pairs if p.similarity >= t), key=lambda p: p.similarity, reverse=True)
