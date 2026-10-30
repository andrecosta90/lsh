import pytest

from lsh.minhash import MinHasher
from lsh.similarity import signature_similarity


def test_signature_has_one_component_per_hash_function():
    assert MinHasher(n=16, seed=1).signature({"ab", "bc"}).shape == (16,)


def test_same_seed_gives_same_signature():
    doc = {"ab", "bc", "cd"}
    assert (MinHasher(n=8, seed=1).signature(doc) == MinHasher(n=8, seed=1).signature(doc)).all()


def test_signature_agreement_estimates_jaccard_similarity():
    a = {f"s{i}" for i in range(0, 100)}
    b = {f"s{i}" for i in range(50, 150)}  # Jaccard = 50 / 150
    minhasher = MinHasher(n=1000, seed=7)

    estimate = signature_similarity(minhasher.signature(a), minhasher.signature(b))

    assert estimate == pytest.approx(1 / 3, abs=0.05)


def test_empty_set_has_no_signature():
    with pytest.raises(ValueError, match="empty"):
        MinHasher(n=4, seed=1).signature(set())
