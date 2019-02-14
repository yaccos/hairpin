import hairpin
import Bio.Seq
import numpy as np

t=Bio.Seq.Seq('ACUGG')

hairpin.find_matches(t, window_size=2)


