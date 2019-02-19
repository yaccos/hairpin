import hairpin
import Bio.Seq
import Bio.SeqIO
import numpy as np

seq = Bio.SeqIO.read("data/sequence.fasta", format='fasta')

t=Bio.Seq.Seq('ACUGG')

hairpin.find_matches(t, window_size=2)


