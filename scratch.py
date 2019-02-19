import Bio.SeqIO
import Bio.Seq
object = Bio.SeqIO.parse('/Users/Jakob/Documents/OneDrive - NTNU/MOL3022/Project/hairpin/data/ls_orchid.fasta','fasta')

for seq in object:
    print(seq.id)
    print(repr(seq.seq))
    print(len(seq))

print(__name__)

import numpy as np
import matplotlib.pyplot as plt
x = np.arange(0.1, 4, 0.5)
y = np.exp(-x)
x_1 = np.array([0, 1, 2])
y_1 = np.array([3, 5, 10])
yerr = 0.1 + 0.2*np.sqrt(x)
xerr = 0.1 + yerr
plt.figure()
plt.errorbar(x, y, xerr=0.2, yerr=0.4,fmt = 'o', label='Randomized results with error bars')
plt.plot(x_1, y_1, 'ro',label='Significant intervals')
plt.title("Simplest errorbars, 0.2 in x, 0.4 in y")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('data/plot.pdf', format='pdf')
plt.show()



A = np.array([1, 2, 3])
B = np.array([3, 4, 5])
np.subtract.outer(A, B)

mymatrix = np.array([[11, 12, 13],
                      [21, 22, 23],
                      [31, 32, 33]])

a,b = mymatrix[[0, 1],]

mymatrix[((1, 2), (1,1))]

def myfunction( x ):
    return sum(x)

temp= myfunction(mymatrix)
np.random.choice([b'A', b'C', b'G', b'U'], size=100, p=[0.2, 0.7, 0.05, 0.05])

transcript = Bio.Seq.Seq('AGUAC')
shifted_transcript = transcript[2:]

np.array(transcript, dtype=np.string_)[:-2] == np.array(shifted_transcript, dtype=np.string_)

window_size=2
transcript = Bio.Seq.Seq('AUCCAG')
matches_in_frame = np.zeros(shape=len(transcript), dtype=np.bool_)
np_transcript = np.array(transcript, dtype=np.string_)
np_shifted_transcript =np_transcript[window_size:]
test_transcript = np_transcript[:-window_size]
a=np.logical_and(test_transcript ==
               b'G',
               np.isin(np_shifted_transcript, (b'C', b'U')))
b=np.logical_and(test_transcript ==
               b'U',
               np.isin(np_shifted_transcript,
                       (b'G', b'A')))
c=np.logical_and(test_transcript ==
               b'A',
               np_shifted_transcript == b'U')
d=np.logical_and(test_transcript == b'C',
               np_shifted_transcript == b'G')
print(a, b, c, d)
np.logical_or.reduce(a,b,c,d)
matches_in_frame[:-window_size] = (test_transcript == b'G') * \
                                  np.isin(np_shifted_transcript, (b'U', b'C')) + \
                                  (test_transcript == b'U') * \
                                  np.isin(np_shifted_transcript, (b'A', b'G')) + \
                                  (test_transcript == b'C') * \
                                  (np_shifted_transcript == b'G') + \
                                  (test_transcript == b'A') * \
                                  (np_shifted_transcript == b'U')