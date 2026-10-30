"""Locality-Sensitive Hashing for documents (MMDS, section 3.4.3)."""

from lsh.banding import candidate_probability, choose_bands, threshold
from lsh.index import LSHIndex
from lsh.minhash import MinHasher
from lsh.pipeline import SimilarPair, find_similar
from lsh.shingling import shingles
from lsh.similarity import jaccard, signature_similarity

__all__ = [
    "LSHIndex",
    "MinHasher",
    "SimilarPair",
    "candidate_probability",
    "choose_bands",
    "find_similar",
    "jaccard",
    "shingles",
    "signature_similarity",
    "threshold",
]
