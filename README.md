Developed by Andre H. Costa Silva
11/09/2013

========================

Source code based on: http://infolab.stanford.edu/~ullman/mmds.html
3.4.3 Combining the Techiniques
Mining of Massive Datasets
Rajaraman,Leskovec,Ullman

========================

1. Pick a value of k and construct from each document the set of k-shingles. Optionally, hash the k-shingles to shorter bucket numbers;for this set 'h=True'. Default h=False.

>> docs,shingles = construct_set_shingles(docs,k)
or >> docs,shingles = construct_set_shingles(docs,k,h=True)

========================
2. Sort the document-shingle pairs to order them by shingle.

>> matrix = sort_document_shingle(docsChange,shingles)

========================
3. Pick a length 'n' for the minHash signatures. That is, pick 'n'randomly chosen hash functions. Default n = 12

>> SIG = compute_minHash_signatures(matrix,n)

========================

4. Choose a threshold t that defines how similar documents have to be.
5. Construct candidate pairs by applying the LSH technique of Section 3.4.1.
6. Examine each candidate pair's signatures and determine whether the fraction of components in which they agree is at least t. Pick a number of bands b and a number of rows r such that br = n, and the threshold t is approximately (1/b)^(1/r). Default b=4 and r=3.

>> candidates,sort = apply_LSH_technique(SIG,0.7,b,r)

========================
