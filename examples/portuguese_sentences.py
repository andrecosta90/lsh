"""The example from the 2013 version, run through the corrected pipeline."""

from lsh import find_similar

DOCS = [
    "eu conheco voce",
    "penso que te conheco",
    "eu fiz aquilo",
    "fiz isso",
    "ela disse que conhece voce",
    "conheco voce pessoalmente",
    "acho que eu conheco voce",
    "eu conheco voce pessoalmente",
]

for pair in find_similar(dict(enumerate(DOCS)), k=3, n=200, t=0.5, seed=0):
    print(f"{pair.similarity:6.1%}  {DOCS[pair.a]!r} ~ {DOCS[pair.b]!r}")
