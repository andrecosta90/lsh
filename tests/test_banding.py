import pytest

from lsh.banding import candidate_probability, choose_bands, threshold


@pytest.mark.parametrize(
    "s, expected",
    # MMDS Figure 3.11: b = 20, r = 5
    [(0.2, 0.006), (0.3, 0.047), (0.4, 0.186), (0.5, 0.470), (0.6, 0.802), (0.7, 0.975), (0.8, 0.9996)],
)
def test_candidate_probability_matches_book_figure_3_11(s, expected):
    assert candidate_probability(s, bands=20, rows=5) == pytest.approx(expected, abs=5e-4)


def test_threshold_of_20_bands_of_5_rows_is_about_0_55():
    # MMDS 3.4.2: "the threshold ... is approximately (1/20)^(1/5) = 0.55"
    assert threshold(bands=20, rows=5) == pytest.approx(0.55, abs=0.005)


@pytest.mark.parametrize("t, expected", [(0.55, (20, 5)), (0.8, (10, 10)), (0.15, (50, 2))])
def test_choose_bands_picks_the_factorisation_of_n_closest_to_t(t, expected):
    # thresholds for n = 100: (50, 2) ≈ .14, (20, 5) ≈ .55, (10, 10) ≈ .79
    assert choose_bands(n=100, t=t) == expected


@pytest.mark.parametrize(
    "n, expected",
    # P(candidate | s = 0.5) for n = 100: (25, 4) ≈ .80, (50, 2) ≈ 1.0
    # and for n = 200: (40, 5) ≈ .72, (50, 4) ≈ .96
    [(100, (50, 2)), (200, (50, 4))],
)
def test_min_recall_picks_the_highest_threshold_that_still_reaches_it(n, expected):
    assert choose_bands(n=n, t=0.5, min_recall=0.95) == expected


def test_unreachable_min_recall_is_an_error():
    # n = 1 only allows b = r = 1, where P(candidate | s = t) = t.
    with pytest.raises(ValueError, match="min_recall"):
        choose_bands(n=1, t=0.5, min_recall=0.95)
