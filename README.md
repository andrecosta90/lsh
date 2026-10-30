# lsh

Finding similar documents with **Locality-Sensitive Hashing**, following section
3.4.3 *Combining the Techniques* of
[Mining of Massive Datasets](http://www.mmds.org/) (Leskovec, Rajaraman, Ullman).

A rewrite of a 2013 attempt, kept in [`legacy/`](legacy/) for comparison.

## Usage

```python
from lsh import find_similar

docs = {"a": "eu conheco voce pessoalmente", "b": "conheco voce pessoalmente", "c": "fiz isso"}
find_similar(docs, k=3, n=200, t=0.8, seed=0)
# [SimilarPair(a='a', b='b', similarity=0.8846...)]
```

```sh
uv run python examples/portuguese_sentences.py
uv run pytest
```

## The seven steps of 3.4.3 and where they live

| Step | Book | Code |
|---|---|---|
| 1. Pick `k` and build each document's set of k-shingles, optionally hashed to 4 bytes | 3.2 | `shingling.shingles`, `minhash._row_number` |
| 2. Sort document–shingle pairs by shingle | 3.3.5 | implicit, see below |
| 3. Pick `n` and compute minhash signatures | 3.3.5 | `MinHasher` |
| 4. Choose threshold `t`, then `b`, `r` with `b·r = n` and `t ≈ (1/b)^(1/r)` | 3.4.2 | `banding.choose_bands`, `threshold`, `candidate_probability` |
| 5. Build candidate pairs with the banding technique | 3.4.1 | `LSHIndex` |
| 6. Keep candidates whose signatures agree in at least a fraction `t` of components | 3.3.3 | `similarity.signature_similarity` |
| 7. Optionally check the documents' real (Jaccard) similarity | 3.1.1 | `similarity.jaccard`, `find_similar(verify=True)` |

`find_similar` chains them all.

**Choosing `b` and `r`:** by default `choose_bands` picks the threshold closest to `t`,
so a pair with similarity close to `t` is missed fairly often (with `n = 200, t = 0.5`
it becomes a candidate only ~72% of the time). When false negatives are what you need
to avoid, pass `min_recall=0.95`: this picks the highest threshold at which a pair of
similarity `t` still becomes a candidate with probability ≥ 0.95, which costs more
false positives for steps 6–7 to filter out. Step 6 itself compares a *noisy estimate*
against `t`, so pairs right at `t` are still dropped about half the time.

**Steps 2 and 3:** the algorithm in 3.3.5 scans the sorted rows of the characteristic
matrix and keeps, for each hash function `h(x) = (a·x + b) mod p`, the smallest
`h(row)` among the rows where the document has a 1. That minimum depends only on the
document's set of shingles, so `MinHasher` computes it straight from the set and
never builds the matrix.

## What was wrong with the 2013 version

- **Shingle hash** ([#1](https://github.com/andrecosta90/lsh/issues/1)): `sum(ord(c)) % 101` sends anagrams to the same value and allows only 101 values.
- **Minhash functions:** used `mod n_rows` with `a ∈ [0, 1000]` (`a` could be 0), so they weren't permutations. Now `mod p` with `p = 2³¹ − 1` prime and `a ≥ 1`.
- **Band hash:** `sum(band) % 101` ignored order and collided often. Bands now have to be identical.
- **Buckets:** only the first two documents in a bucket were paired. Now every pair is.
- **Step 6:** compared signatures with cosine similarity. The book's measure is the fraction of agreeing components.
- **Steps 4 and 7:** `b` and `r` weren't derived from `t`, and the documents were never checked.
