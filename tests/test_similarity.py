import numpy as np
import pytest

from lsh.similarity import jaccard, signature_similarity


def test_jaccard_is_intersection_over_union():
    assert jaccard({"a", "b", "c", "d"}, {"c", "d", "e"}) == pytest.approx(2 / 5)


def test_signature_similarity_is_fraction_of_agreeing_components():
    assert signature_similarity(np.array([1, 2, 3, 4]), np.array([1, 0, 3, 0])) == 0.5
