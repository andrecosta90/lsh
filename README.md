

# LSH - Locality Sensitive Hashing Implementation

## Overview
A Python implementation of **Locality Sensitive Hashing (LSH)** for finding similar documents. This script follows **Section 3.4.3 ("Combining the Techniques")** of the textbook *Mining of Massive Datasets* by Leskovec, Rajaraman, and Ullman. It processes raw text documents through shingling, MinHash signature generation, and band-based LSH to efficiently identify document pairs above a specified cosine similarity threshold.

## Features
- **k-Shingle Extraction**: Splits documents into fixed-length character sequences (`k`)
- **Character Matrix Construction**: Builds a binary document-shingle matrix
- **MinHash Signatures**: Computes compact signatures using randomized hash functions
- **LSH Banding**: Partitions signatures into bands and rows to generate candidate pairs
- **Similarity Filtering**: Validates candidates using cosine distance against a threshold `t`

## Prerequisites
- Python 2.7 *(Note: The code uses Python 2 syntax such as `xrange`, classic `print` statements, and `raise 'string'`)*
- `numpy`

## Installation
Install the required dependency:
```bash
pip install numpy
```

## Usage
Run the script directly from the terminal:
```bash
python LSH_impl.py
```
The script contains a built-in example using a list of Portuguese sentences. Upon execution, it will construct shingles, generate MinHash signatures, apply LSH, and print the top `j` most similar document pairs along with their similarity percentages.

## Configuration
Algorithm parameters can be adjusted at the bottom of `LSH_impl.py`:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `k`       | Shingle character length | `3` |
| `n`       | Number of hash functions (signature length) | `1000` |
| `t`       | Cosine similarity threshold | `0.8` |
| `bands`   | Number of LSH bands | `100` |
| `rows`    | Rows per band | `10` |

⚠️ **Important**: LSH parameters must satisfy `bands * rows == n`. The approximate threshold is `(1/bands)^(1/rows)`.

## Python 3 Compatibility
If running on Python 3, you will need to:
1. Replace `xrange` with `range`
2. Update `print` statements to use function syntax: `print(...)`
3. Change `raise 'error message'` to `raise Exception('error message')`

## References
- Leskovec, J., Rajaraman, A., & Ullman, J. D. *Mining of Massive Datasets*, Chapter 3.4.3
- Original implementation by Andre H. Costa Silva (2013)
