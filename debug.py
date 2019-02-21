import hairpin
import Bio.SeqIO
import numpy as np
infile = 'data/e_coli_short.fasta'
min_window = 4
n_randomizations = 100
max_window = 50
min_length = 25
max_length = 200
q_crit = 0.05
min_motif_length = 7
chrom = 'chr'
chromstart = 0
res_prefix = ''
np.random.seed(334)
sequence = Bio.SeqIO.read(infile, 'fasta', alphabet=Bio.Alphabet.generic_dna).seq
results = hairpin.predict_RNA_genes(sequence, n_randomizations=n_randomizations, min_window=min_window,
                                 max_window=max_window, min_length=min_length,
                                 max_length=max_length, q_crit=q_crit, min_motif_length=min_motif_length)
hairpin.write_results(results=results, chrom=chrom, chromStart=chromstart, res_prefix=res_prefix)
