#!/usr/bin/env python3
import time
import hairpin
import Bio.Seq
import Bio.SeqIO
import Bio.Alphabet
import numpy as np
infile = 'data/e_coli_k_12.fasta'
min_window = 4
n_randomizations = 2000
max_window = 50
min_length = 25
max_length = 200
q_crit = 0.01
min_motif_length = 8
chrom = 'chr'
chromstart = 0
res_prefix = 'shuffle_run_2'
np.random.seed(684640245)
start = time.time()
sequence = Bio.SeqIO.read(infile, 'fasta', alphabet=Bio.Alphabet.generic_dna).seq
shuffeled_sequence_list = list(sequence)
np.random.shuffle(shuffeled_sequence_list)
shuffeled_sequence = Bio.Seq.Seq(data=''.join(shuffeled_sequence_list), alphabet=Bio.Alphabet.generic_dna)
results = hairpin.predict_RNA_genes(shuffeled_sequence, n_randomizations=n_randomizations, min_window=min_window,
                                 max_window=max_window, min_length=min_length,
                                 max_length=max_length, q_crit=q_crit, min_motif_length=min_motif_length)
hairpin.write_results(results=results, chrom=chrom, chromStart=chromstart, res_prefix=res_prefix)
end = time.time()
print('Elapsed time:', end - start)
