import pytest

from lsh.pipeline import SimilarPair, find_similar

# The example from the 2013 version.
DOCS = {
    0: "eu conheco voce",
    1: "penso que te conheco",
    2: "eu fiz aquilo",
    3: "fiz isso",
    4: "ela disse que conhece voce",
    5: "conheco voce pessoalmente",
    6: "acho que eu conheco voce",
    7: "eu conheco voce pessoalmente",
}


def test_finds_the_only_pair_above_the_threshold():
    # True Jaccard of 3-shingles: (5, 7) = 23/26 ≈ 0.885; the next best pair is 0.591.
    assert find_similar(DOCS, k=3, n=200, t=0.8, seed=0) == [SimilarPair(5, 7, pytest.approx(23 / 26))]


def test_without_verification_similarity_is_estimated_from_signatures():
    [pair] = find_similar(DOCS, k=3, n=200, t=0.8, seed=0, verify=False)

    assert (pair.a, pair.b) == (5, 7)
    assert (pair.similarity * 200).is_integer()  # a fraction of the 200 components
    assert pair.similarity == pytest.approx(23 / 26, abs=0.1)


def test_lowercase_makes_documents_differing_only_in_case_identical():
    docs = {"a": "Eu conheco voce", "b": "eu conheco voce"}

    assert find_similar(docs, k=3, n=100, t=0.9, seed=0) == []
    assert find_similar(docs, k=3, n=100, t=0.9, seed=0, lowercase=True) == [SimilarPair("a", "b", 1.0)]


def test_results_are_ordered_from_most_to_least_similar():
    docs = {"x": "abcdefghij", "y": "abcdefghij", "z": "abcdefghiX"}  # x, y vs z: 7/9

    result = find_similar(docs, k=3, n=200, t=0.6, seed=0)

    assert [p.similarity for p in result] == pytest.approx([1.0, 7 / 9, 7 / 9])


def test_min_recall_is_used_to_choose_bands():
    with pytest.raises(ValueError, match="min_recall"):
        find_similar(DOCS, k=3, n=1, t=0.5, seed=0, min_recall=0.95)


def test_documents_shorter_than_k_are_reported_by_key():
    docs = {"ok": "eu conheco voce", "d3": "oi", "d9": ""}

    with pytest.raises(ValueError, match=r"shorter than k=3: \['d3', 'd9'\]"):
        find_similar(docs, k=3, n=100, t=0.5, seed=0)
