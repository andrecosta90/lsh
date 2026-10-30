"""Step 1 of MMDS 3.4.3: represent each document as its set of k-shingles (3.2)."""


def shingles(text: str, k: int, *, lowercase: bool = False) -> frozenset[str]:
    """Return the set of all substrings of length ``k`` found in ``text``.

    Any run of whitespace is treated as a single blank (MMDS 3.2.2).
    """
    if k < 1:
        raise ValueError(f"k must be positive, got {k}")
    text = " ".join(text.split())
    if lowercase:
        text = text.lower()
    return frozenset(text[i : i + k] for i in range(len(text) - k + 1))
