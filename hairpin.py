#!/usr/bin/env python3
import sys
import Bio
import Bio.SeqRecord
import Bio.SeqIO
import Bio.Alphabet
import numpy as np
import scipy.stats
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


def count_consecutive_matches(match_array, min_motif_length=4):
    """
    :param match_array: An boolean Numpy array, showing matches or non-matches
    :param min_motif_length: The minimum length of a match (consecutive trues in the array) in order to report it
    :return: Two integer Numpy arrays of the same length,
    res_start (first return element) reports the position in which a significant motif starts,
    while res_length (second return element) reports the corresponding length of the motif
    """
    # Pads the sequence with False at each end in order to detect consecutive occurrences of True values
    # (correspond to nucleotide matches)
    padded_array = np.concatenate([False], match_array, [False])
    false_indecies = np.nonzero(np.invert(padded_array))
    n_consecutive_trues = np.diff(false_indecies) - 1
    # Finds the consecutive sequences of base pairings exceeding our threshold
    significant_trues = n_consecutive_trues >= min_motif_length
    # Number of significant motifs
    n_significant_sequences = np.sum(significant_trues)
    # Length of the moftis
    res_length = n_consecutive_trues[significant_trues]
    # The start positions of the significant motifs
    res_start = false_indecies[significant_trues]
    return res_start, res_length


class MatchSequence:
    def __init__(self, match_array, start_nucleotide, start_window, min_motif_length=4):
        self.match_array = match_array
        self.start_nucleotide = start_nucleotide
        self.start_window = start_window
        self.min_motif_length = min_motif_length

    def find_significant_motifs(self):
        res_start, res_length = count_consecutive_matches(self.match_array, min_motif_length=self.min_motif_length)
        return self._translate(res_start=res_start, res_length=res_length)

    def translate(self, res_start, res_length):
        """
     Given a result from count_consecutive_motifs, this function relates the significant hitfs onto the original
     sequence
     :param res_start:
     :param res_length:
     :return:
     """
        start_position = self.start_nucleotide + res_start
        match_length = res_length
        end_position = self.start_nucleotide + self.start_window - 2 * self.start_nucleotide
        return start_position, match_length, end_position


def motif_bases(start_position, motif_length, end_position):
    first_part = np.array(range(start_position, start_position + motif_length))
    secound_part = np.array(range(end_position, end_position - motif_length), -1)
    return np.concatenate([first_part, secound_part])


class MatchResult:
    def __init__(self, start_position, motif_length, end_position):
        self.start_position = start_position
        self.motif_length = motif_length
        self.end_position = end_position
        self.order_start_position = np.argsort(start_position)
        self.order_end_position = np.argsort(end_position)

    def __len__(self):
        return len(self.start_position)

    def count_matched_bases(self, start_motif, end_motif):
        """
        :param start_motif: Index of the motif (relative to the order of start_position),
         marking the beginning the region to be searched
        :param end_motif: Index of the motif (relative to the order of end_position),
         marking the end of the region to be search
        :return: Number of matched bases in the region, counting only the motifs being entirely inside the region
        """
        start = self.start_position[self.order_start_position[start_motif]]
        end = self.end_position[self.order_end_position[end_motif]]
        # Finds the motifs contained in the region
        motifs_in_region = np.nonzero((self.start_position >= start) * (self.end_position <= end))
        matched_bases = [motif_bases(self.start_position[motif], self.motif_length[motif], self.end_position[motif])
                         for motif in motifs_in_region]
        matched_bases = np.concatenate(matched_bases)
        count = np.unique(matched_bases).size
        return count


def extract_match_results(parity_match_matrix, max_window):
    n_row = parity_match_matrix.size(axis=0)
    n_col = parity_match_matrix.size(axis=1)
    # Preallocates the arrays storing the results
    start_position = [None for _ in range(n_col + n_row - 1)]
    motif_length = [None for _ in range(n_col + n_row - 1)]
    end_position = [None for _ in range(n_col + n_row - 1)]
    # Extracts the diagonals above and including the main diagonal
    for col in range(n_col):
        start_nucleotide = col
        start_window = max_window
        this_match_diagonal = np.diagonal(parity_match_matrix, col)
        start_position[col], motif_length[col], end_position[col] = \
            MatchSequence(this_match_diagonal, start_nucleotide, start_window).find_significant_motifs()
    # Extracts the diagonals below the main diagonal
    for row in range(1, n_row):
        start_nucleotide = 0
        start_window = max_window - 2 * row
        this_match_diagonal = np.diagonal(parity_match_matrix, -row)
        start_position[n_col + row - 1], motif_length[n_col + row - 1], end_position[n_col + row - 1] = \
            MatchSequence(this_match_diagonal, start_nucleotide, start_window).find_significant_motifs()
    # Collapses the results for all result arrays
    start_position = np.concatenate(start_position)
    motif_length = np.concatenate(motif_length)
    end_position = np.concatenate(end_position)
    return start_position, motif_length, end_position


def find_hairpins(transcript: Bio.SeqRecord.SeqRecord, min_window, max_window, min_motif_length) -> MatchResult:
    # Finds basepair matches for each distance. The first basepair in the match will be reported as TRUE
    matches = [find_matches(transcript, window_size=window_size) for window_size in range(min_window, max_window + 1)]
    # Finds out which window size which are even
    is_even_window = np.array(range(min_window, max_window + 1)) % 2 == 0
    # Constructs a matrix of the matches and reverses the window sizes
    match_matrix = np.vstack(matches)[::-1, :]
    # We separate the odd and even window sizes due to the fact that we want to
    even_match_matrix = match_matrix[np.flipud(is_even_window), :]
    odd_match_matrix = match_matrix[np.invert(np.flipud(is_even_window)), :]
    max_even_window = max_window if max_window % 2 == 0 else max_window - 1
    max_odd_window = max_window if max_window % 2 != 0 else max_window - 1
    even_start_position, even_match_length, even_end_position = extract_match_results(even_match_matrix,
                                                                                      max_even_window)
    odd_start_position, odd_match_length, odd_end_position = extract_match_results(odd_match_matrix,
                                                                                   max_odd_window)
    start_position = np.concatenate([even_start_position, odd_start_position])
    match_length = np.concatenate([even_match_length, odd_match_length])
    end_position = np.concatenate([even_end_position, odd_end_position])
    return MatchResult(start_position, match_length, end_position)


class RandomizationResult:
    def __init__(self, result_matrix: np.ndarray, max_length: int):
        self.mean_matched_bases = np.mean(result_matrix, axis=0)
        self.std_error_matched_bases = np.std(result_matrix, axis=0)
        self.max_length = max_length

    def compare(self, match_result: MatchResult, len_sequence: int, p_crit: int = 0.05):
        # Finds the motif intervals to consider
        length_of_interval = np.subtract.outer(match_result.start_position, match_result.end_position)
        # Finds the feasable intervals
        intervals_to_consider = np.nonzero(np.min(match_result.motif_length) <= length_of_interval <= self.max_length)
        interval_indecies = np.transpose(intervals_to_consider)
        interval_list = map(lambda x: x, interval_indecies)
        n_matched_bases_in_intervals = np.array(map(lambda start, end:
                                                    match_result.count_matched_bases(start_motif=start, end_motif=end),
                                                    interval_list))
        length_of_intervals_to_consider = length_of_interval[intervals_to_consider]
        z_values_for_intervals = (n_matched_bases_in_intervals
                                  - self.mean_matched_bases[length_of_intervals_to_consider]) / \
                                 self.std_error_matched_bases[length_of_intervals_to_consider]
        p_values_for_intervals = 1 - scipy.stats.norm.cdf(z_values_for_intervals)
        correced_p_values_for_intervals = np.power(p_values_for_intervals,
                                                   length_of_intervals_to_consider / len_sequence)
        significant_intervals = np.nonzero(correced_p_values_for_intervals < p_crit)




def compute_randomized_structures(transcript: Bio.SeqRecord.SeqRecord, n_randomizations: int, min_window: int,
                                  max_window: int,
                                  max_length: int, min_motif_length: int):
    raw_weights = np.array([transcript.count('A'), transcript.count('C'), transcript.count('G'), transcript.count('U')])
    weights = raw_weights / np.sum(raw_weights)
    res = [None for _ in range(n_randomizations)]

    def counting_wrapper(transcript, min_window, max_window, min_motif_lenght):
        motif_results = find_hairpins(transcript=transcript, min_window=min_window,
                                      max_window=max_window, min_motif_length=min_motif_length)
        return motif_results.count_matched_bases(start_motif=0, end_motif=len(motif_results))

    for i in range(n_randomizations):
        # Creates a random sequence according to the original distribution in the transcript
        random_sequence = np.random.choice([b'A', b'C', b'G', b'U'], size=max_length, p=weights)
        matched_base_count = [counting_wrapper(transcript=random_sequence[:length], min_window=min_window,
                                               max_window=max_window, min_motif_length=min_motif_length)
                              for length in range(max_length)]
        res[i] = np.array(matched_base_count)
    res_combined = np.vstack(res)
    return RandomizationResult(result_matrix=res_combined, len_sequence = len(transcript), max_length=max_length)


def predict_RNA_genes(sequence, n_randomizations, min_window, max_window, max_length, p_crit, min_motif_length):
    sequence = sequence.upper()
    transcript = sequence.upper()
    randomized_results = compute_randomized_structures(transcript, n_randomizations, min_window,
                                                       max_window, max_length, min_motif_length)
    actual_results = find_hairpins(transcript, min_window, max_window, max_length, min_motif_length)

    return randomized_results.compare(actual_results, p_crit=p_cript)


@click.command()
@click.argument(u'infile')
@click.option('--n_randomizations', default=50, help="Number of randomizations used to assess the significance of the"
                                                     " results")
@click.option('--min_window', default=4, help='Minimal window searching for matching nucleotides')
@click.option('--max_window', default=100, help='Maximal window searching for matching nucleotides')
@click.option('--max_length', default=100, help='Maximal length of reported RNA coding region')
@click.option('--min_motif_length', default=4,
              help='The minimum number of consecutive base parings necessary to take such a motif into '
                   'consideration')
@click.option('--p_crit', default=0.05, help='Critical p-value for the results')
def main(INFILE, min_window, n_randomizations, max_window, max_length, p_crit, min_motif_length):
    """
    :param INFILE: The filename (FASTA format) of the coding DNA strand being used to predict the RNA genes.
    """
    sequences = Bio.SeqIO.parse(infile, 'fasta', alphabet=Bio.Alphabet.generic_dna)
    results = [predict_RNA_genes(sequence, n_randomizations=n_randomizations, min_window=min_window,
                                 max_window=max_window,
                                 max_length=max_length, p_crit=p_crit, min_motif_length=min_motif_length)
               for sequence in sequences]


if __name__ == '__main__':
    if sys.version_info.major == 3 and sys.version_info.minor >= 5:
        main()
    else:
        print('Please use Python 3.5 or later')
