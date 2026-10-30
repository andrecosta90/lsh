import pytest

from lsh.shingling import shingles


def test_k_shingles_of_book_example_3_3():
    assert shingles("abcdabd", k=2) == {"ab", "bc", "cd", "da", "bd"}


def test_any_run_of_whitespace_counts_as_a_single_blank():
    assert shingles("a \t\n b", k=2) == {"a ", " b"}


def test_case_is_preserved_unless_lowercase_is_requested():
    assert shingles("Ab", k=2) == {"Ab"}
    assert shingles("Ab", k=2, lowercase=True) == {"ab"}


def test_text_shorter_than_k_has_no_shingles():
    assert shingles("ab", k=3) == frozenset()


@pytest.mark.parametrize("k", [0, -1])
def test_k_must_be_positive(k):
    with pytest.raises(ValueError):
        shingles("abc", k=k)
