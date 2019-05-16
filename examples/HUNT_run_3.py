#!/usr/bin/env python3
import time
import hairpin
import Bio.SeqIO
import numpy as np
infile = 'data/e_coli_k_12.fasta'
min_window = 4
n_randomizations = 2000
max_window = 50
min_length = 200
max_length = 1000
q_crit = 0.01
min_motif_length = 8
chrom = 'chr'
chromstart = 0
res_prefix = 'k_12_run3'
np.random.seed(684640245)
start = time.time()
sequence = Bio.SeqIO.read(infile, 'fasta', alphabet=Bio.Alphabet.generic_dna).seq
results = hairpin.predict_RNA_genes(sequence, n_randomizations=n_randomizations, min_window=min_window,
                                 max_window=max_window, min_length=min_length,
                                 max_length=max_length, q_crit=q_crit, min_motif_length=min_motif_length)
hairpin.write_results(results=results, chrom=chrom, chromStart=chromstart, res_prefix=res_prefix)
end = time.time()
print('Elapsed time:', end - start)
