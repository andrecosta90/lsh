import numpy as np
import pytest

from lsh.index import LSHIndex


def index_of(bands, rows, **signatures):
    index = LSHIndex(bands=bands, rows=rows)
    for key, sig in signatures.items():
        index.add(key, np.array(sig))
    return index


def test_documents_identical_in_at_least_one_band_are_candidates():
    index = index_of(2, 2, a=[1, 2, 3, 4], b=[1, 2, 9, 9], c=[5, 6, 7, 8])
    assert index.candidate_pairs() == {("a", "b")}


def test_band_order_matters():
    # The 2013 version hashed sum(band), so [1, 2] and [2, 1] collided.
    index = index_of(1, 2, a=[1, 2], b=[2, 1])
    assert index.candidate_pairs() == set()


def test_every_pair_within_a_bucket_is_a_candidate():
    # The 2013 version only paired the first two documents of each bucket.
    index = index_of(1, 2, a=[1, 1], b=[1, 1], c=[1, 1])
    assert index.candidate_pairs() == {("a", "b"), ("a", "c"), ("b", "c")}


def test_signature_length_must_be_bands_times_rows():
    with pytest.raises(ValueError, match="bands \\* rows"):
        index_of(2, 2, a=[1, 2, 3])
