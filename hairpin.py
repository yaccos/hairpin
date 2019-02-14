import Bio
import Bio.SeqRecord
import Bio.SeqIO
import Bio.Alphabet
import numpy as np
import click


def find_matches(transcript: Bio.SeqRecord.SeqRecord, window_size) -> np.array:
    """
    :param transcript: The RNA transcipt of be analyzed
    :param window_size: The window size for which the matching is done
    :return: A boolean Numpy-array where the first residue participating in a bound at distance corresponding to the
    window size is marked as True
    Note that when a potential bound is found, only the first residue of that pair is marked
    """

    matches_in_frame = np.zeros(shape=len(transcript), dtype=np.bool_)
    np_transcript = np.array(list(transcript), dtype=np.string_)
    np_shifted_transcript = np_transcript[window_size:]
    test_transcript = np_transcript[:-window_size]

    # To account for the consensus base pairs {G,C}, {A,U} and {G,U} in a vectorized manner
    matches_in_frame[:-window_size] = (test_transcript == b'G') * \
                                      np.isin(np_shifted_transcript, (b'U', b'C')) + \
                                      (test_transcript == b'U') * \
                                      np.isin(np_shifted_transcript, (b'A', b'G')) + \
                                      (test_transcript == b'C') * \
                                      (np_shifted_transcript == b'G') + \
                                      (test_transcript == b'A') * \
                                      (np_shifted_transcript == b'U')

    return matches_in_frame


def count_consecutive_matches(match_array, min_match_length=4):
    """
    :param match_array: An boolean Numpy array, showing matches or non-matches
    :param min_match_length: The minimum length of a match (consecutive trues in the array) in order to report it
    :return: Two integer Numpy arrays of the same length,
    res_start (first return element) reports the position in which a significant matching starts,
    while res_length (second return element) reports the corresponding length of the matchings
    """
    # Pads the sequence with False at each end in order to detect consecutive occurrences of True values
    # (correspond to nucleotide matches)
    padded_array = np.concatenate([False], match_array, [False])
    false_indecies = np.nonzero(np.invert(padded_array))
    n_consecutive_trues = np.diff(false_indecies) - 1
    # Finds the consecutive sequences of base pairings exceeding our threshold
    significant_trues = n_consecutive_trues >= min_match_length
    # Number of significant matches
    n_significant_sequences = np.sum(significant_trues)
    # Length of the matchings
    res_length = n_consecutive_trues[significant_trues]
    # The start positions of the significant matchings
    res_start = false_indecies[significant_trues]
    return res_start, res_length


class MatchSequence:
    def __init__(self, match_array, start_nucleotide, start_window, min_match_length=4):
        self.match_array = match_array
        self.start_nucleotide = start_nucleotide
        self.start_window = start_window
        self.min_match_length = min_match_length

    def find_significant_matches(self):
        res_start, res_length = count_consecutive_matches(self.match_array, min_match_length=self.min_match_length)
        return self._translate(res_start=res_start,res_length=res_length)


    def translate(self, res_start, res_length):
        """
     Given a result from count_consecutive_matches, this function relates the significant hits onto the original
     sequence
     :param res_start:
     :param res_length:
     :return:
     """
        start_position = self.start_nucleotide + res_start
        match_length = res_length
        end_position = self.start_nucleotide + self.start_window - 2 * self.start_nucleotide
        return (start_position, match_length, end_position)


def extract_match_results(parity_match_matrix, max_window):
    n_row = parity_match_matrix.size(axis=0)
    n_col = parity_match_matrix.size(axis=1)
    res = []
    for col in range(n_col):
        start_nucleotide = col
        start_window = max_window

        


def find_hairpins(sequence: Bio.SeqRecord.SeqRecord, min_window, max_window, min_match_length):
    transcript = sequence.transcribe()
    # Finds basepair matches for each distance. The first basepair in the match will be reported as TRUE
    matches = [find_matches(transcript, window_size=window_size) for window_size in range(min_window, max_window + 1)]
    # Finds out which window size which are even
    is_even_window = np.array(range(min_window, max_window + 1)) % 2 == 0
    # Constructs a matrix of the matches and reverses the window sizes
    match_matrix = np.vstack(matches)[::-1, :]
    # We separate the odd and even window sizes due to the fact that we want to
    even_match_matrix = match_matrix[np.flipud(is_even_window), :]
    odd_match_matrix = match_matrix[np.invert(np.flipud(is_even_window)), :]
    for i in even_match_matrix.size(axis=0):
        match_sequence = MatchSequence()


def predict_RNA_genes(sequence, max_window, max_length, q_crit):
    sequence = sequence.upper()


@click.command()
@click.argument(u'infile')
@click.option('--min_window', default=4, help='Minimal window searching for matching nucleotides')
@click.option('--max_window', default=100, help='Maximal window searching for matching nucleotides')
@click.option('--max_length', default=100, help='Maximal length of reported RNA coding region')
@click.option('--q_crit', default=0.05, help='Critical q-value for the results')
def main(infile, min_window, max_window, max_length, q_crit):
    """
    :param INFILE: The filename (FASTA format) of the coding DNA strand being used to predict the RNA genes.
    """

    sequences = Bio.SeqIO.parse(infile, 'fasta', alphabet=Bio.Alphabet.generic_dna)
    results = [predict_RNA_genes(sequence, max_window, max_length, q_crit) for sequence in sequences]


if __name__ == '__main__':
    main()
