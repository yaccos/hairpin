#!/usr/bin/env python3
import sys
import Bio
import Bio.SeqRecord
import Bio.SeqIO
import Bio.Seq
import Bio.Alphabet
import numpy as np
import scipy.stats
import click
import joblib
import matplotlib.pyplot as plt
import statsmodels.stats.multitest


def find_matches(np_transcript: np.ndarray, window_size: int) -> np.ndarray:
    """
    :param np_transcript: The RNA transcipt of be analyzed
    :param window_size: The window size for which the matching is done
    :return: A boolean Numpy-array where the first residue participating in a bound at distance corresponding to the
    window size is marked as True
    Note that when a potential bound is found, only the first residue of that pair is marked
    """

    matches_in_frame = np.zeros(shape=np_transcript.size, dtype=np.bool_)
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


def count_consecutive_matches(match_array: np.ndarray, min_motif_length=4):
    """
    :param match_array: An boolean Numpy array, showing matches or non-matches
    :param min_motif_length: The minimum length of a match (consecutive trues in the array) in order to report it
    :return: Two integer Numpy arrays of the same length,
    res_start (first return element) reports the position in which a significant motif starts,
    while res_length (second return element) reports the corresponding length of the motif
    """
    # Pads the sequence with False at each end in order to detect consecutive occurrences of True values
    # (correspond to nucleotide matches)
    if len(match_array) == 0:
        return np.array([]), np.array([])
    padded_array = np.concatenate([[False], match_array, [False]])
    false_indecies = np.nonzero(np.invert(padded_array))[0]
    n_consecutive_trues = np.diff(false_indecies) - 1
    # Finds the consecutive sequences of base pairings exceeding our threshold
    significant_trues = n_consecutive_trues >= min_motif_length
    # Number of significant motifs
    n_significant_sequences = np.sum(significant_trues)
    # Length of the moftis
    res_length = n_consecutive_trues[significant_trues]
    # The start positions of the significant motifs
    res_start = false_indecies[:-1][significant_trues]
    return res_start, res_length


class MatchSequence:
    def __init__(self, match_array, start_nucleotide, start_window, min_motif_length=4):
        self.match_array = match_array
        self.start_nucleotide = start_nucleotide
        self.start_window = start_window
        self.min_motif_length = min_motif_length

    def find_significant_motifs(self):
        res_start, res_length = count_consecutive_matches(self.match_array, min_motif_length=self.min_motif_length)
        return self.translate(res_start=res_start, res_length=res_length)

    def translate(self, res_start: np.ndarray, res_length: np.ndarray):
        """
     Given a result from count_consecutive_motifs, this function relates the significant motifs onto the original
     sequence
     :param res_start:
     :param res_length:
     :return:
     """
        start_position = self.start_nucleotide + res_start
        match_length = res_length
        end_position = self.start_nucleotide + res_start + self.start_window - 2 * res_start
        return start_position, match_length, end_position


def motif_bases(start_position, motif_length, end_position):
    first_part = np.array(range(start_position, start_position + motif_length))
    secound_part = np.array(range(end_position, end_position - motif_length, -1))
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

    def count_matched_bases(self, start_motif, end_motif) -> int:
        """
        :param start_motif: Index of the motif,
         marking the beginning the region to be searched
        :param end_motif: Index of the motif,
         marking the end of the region to be search
        :return: Number of matched bases in the region, counting only the motifs being entirely inside the region
        """
        if len(self) == 0:
            return 0
        start = self.start_position[start_motif]
        end = self.end_position[end_motif]
        # Finds the motifs contained in the region
        motifs_in_region = np.nonzero((self.start_position >= start) * (self.end_position <= end))[0]
        matched_bases = [motif_bases(self.start_position[motif], self.motif_length[motif], self.end_position[motif])
                         for motif in motifs_in_region]
        if len(matched_bases) == 0:
            return 0
        matched_bases = np.concatenate(matched_bases)
        count = np.unique(matched_bases).size
        return count


def extract_match_results(parity_match_matrix, max_window, min_motif_length=4):
    n_row = parity_match_matrix.shape[0]
    n_col = parity_match_matrix.shape[1]
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
            MatchSequence(this_match_diagonal, start_nucleotide,
                          start_window, min_motif_length=min_motif_length).find_significant_motifs()
    # Extracts the diagonals below the main diagonal
    for row in range(1, n_row):
        start_nucleotide = 0
        start_window = max_window - 2 * row
        this_match_diagonal = np.diagonal(parity_match_matrix, -row)
        start_position[n_col + row - 1], motif_length[n_col + row - 1], end_position[n_col + row - 1] = \
            MatchSequence(this_match_diagonal, start_nucleotide,
                          start_window, min_motif_length=min_motif_length).find_significant_motifs()
    # Collapses the results for all result arrays
    start_position = np.concatenate(start_position)
    motif_length = np.concatenate(motif_length)
    end_position = np.concatenate(end_position)
    return start_position, motif_length, end_position


def find_hairpins(np_transcript: np.ndarray, min_window: int, max_window: int, min_motif_length: int,
                  parallel: bool = True) -> MatchResult:
    # Finds basepair matches for each distance. The first basepair in the match will be reported as TRUE
    if parallel:
        matches = joblib.Parallel(n_jobs=joblib.cpu_count())(joblib.delayed(find_matches)
                                                             (np_transcript, window_size=window_size)
                                                             for window_size in range(min_window, max_window + 1))
    else:
        matches = [find_matches(np_transcript, window_size=window_size)
                   for window_size in range(min_window, max_window + 1)]
    # Finds out which window size which are even
    is_even_window = np.array(range(min_window, max_window + 1)) % 2 == 0
    # Constructs a matrix of the matches and reverses the window sizes
    match_matrix = np.vstack(matches)[::-1, :]
    # We separate the odd and even window sizes due to the fact that we want to
    even_match_matrix = match_matrix[np.flipud(is_even_window), :]
    odd_match_matrix = match_matrix[np.invert(np.flipud(is_even_window)), :]
    max_even_window = max_window if max_window % 2 == 0 else max_window - 1
    max_odd_window = max_window if max_window % 2 != 0 else max_window - 1
    if parallel:
        (even_start_position, even_match_length, even_end_position), \
        (odd_start_position, odd_match_length, odd_end_position) = joblib.Parallel(n_jobs=2)(
            joblib.delayed(extract_match_results)(call[0], call[1], min_motif_length=min_motif_length)
            for call in ((even_match_matrix, max_even_window), (odd_match_matrix, max_odd_window))
        )
    else:
        even_start_position, even_match_length, even_end_position = \
            extract_match_results(even_match_matrix,
                                  max_even_window,
                                  min_motif_length=min_motif_length)
        odd_start_position, odd_match_length, odd_end_position = \
            extract_match_results(odd_match_matrix,
                                  max_odd_window,
                                  min_motif_length=min_motif_length)
    start_position = np.concatenate([even_start_position, odd_start_position])
    match_length = np.concatenate([even_match_length, odd_match_length])
    end_position = np.concatenate([even_end_position, odd_end_position])
    return MatchResult(start_position, match_length, end_position)


class RandomizationResult:
    def __init__(self, result_matrix: np.ndarray, min_length: int, max_length: int):
        self.mean_matched_bases = np.mean(result_matrix, axis=0)
        self.std_error_matched_bases = np.std(result_matrix, axis=0)
        self.min_length = min_length
        self.max_length = max_length

    def __len__(self):
        return self.mean_matched_bases.size

    def compare(self, match_result: MatchResult, len_sequence: int, min_window: int,
                min_motif_length: int, q_crit: float = 0.05):
        # Finds the motif intervals to consider
        length_of_interval = np.subtract.outer(match_result.end_position, match_result.start_position) + 1
        intervals_to_consider = np.nonzero(np.logical_and(self.min_length <= length_of_interval,
                                                          length_of_interval <= self.max_length))
        interval_indices = np.transpose(intervals_to_consider)
        interval_list = map(lambda x: x, interval_indices)
        n_matched_bases_in_intervals_list = joblib.Parallel(n_jobs=joblib.cpu_count())(joblib.delayed(lambda interval:
                                                                                                      match_result.count_matched_bases(
                                                                                                          start_motif=
                                                                                                          interval[1],
                                                                                                          end_motif=
                                                                                                          interval[0]))(
            interval)
                                                                                       for interval in interval_list)
        n_matched_bases_in_intervals = np.array(n_matched_bases_in_intervals_list)
        length_of_intervals_to_consider = length_of_interval[intervals_to_consider]
        z_values_for_intervals = (n_matched_bases_in_intervals
                                  - self.mean_matched_bases[length_of_intervals_to_consider - self.min_length]) / \
                                 self.std_error_matched_bases[length_of_intervals_to_consider - self.min_length]
        # Discards the intervals where the randomizations do not report any variation in the number of base pairings
        z_values_for_intervals[self.std_error_matched_bases[length_of_intervals_to_consider -
                                                            self.min_length] < np.finfo(float).eps] = -np.Inf
        p_values = 1 - scipy.stats.norm.cdf(z_values_for_intervals)
        p_values = 1 - np.power(1 - p_values, len_sequence / length_of_intervals_to_consider)
        end_motif, start_motif = interval_indices[:, 0], interval_indices[:, 1]
        start_nucleotide = match_result.start_position[start_motif]
        end_nucleotide = match_result.end_position[end_motif]
        p_values_order = np.argsort(p_values)
        res = IntervalResult(randomized_result=self)
        # Removes the overlaps, starting with the most significant interval
        # Indecies relative to the array significant_intervals
        discard = np.zeros(shape=p_values.size, dtype=np.bool_)
        for i in range(p_values.size):
            this_index = p_values_order[i]
            if discard[this_index]:
                continue
            # Adds the interval to the results
            res.add(start_nucleotide=start_nucleotide[this_index], end_nucleotide=end_nucleotide[this_index],
                    p_value=p_values[this_index],
                    n_matched_bases=n_matched_bases_in_intervals[this_index])
            discard[(end_nucleotide >= start_nucleotide[this_index]) *
                    (start_nucleotide <= end_nucleotide[this_index])] = True
        res.correct(q_crit=q_crit)
        return res


class IntervalResult:
    def __init__(self, randomized_result: RandomizationResult):
        self.start_nucleotide = []
        self.end_nucleotide = []
        self.p_value = []
        self.n_matched_bases = []
        self.randomized_result = randomized_result

    def add(self, start_nucleotide, end_nucleotide, p_value, n_matched_bases):
        self.start_nucleotide.append(start_nucleotide)
        self.end_nucleotide.append(end_nucleotide)
        self.p_value.append(p_value)
        self.n_matched_bases.append(n_matched_bases)

    def correct(self, q_crit):
        # Benjamini-Yekutieli correction of p-values
        significant, q_value, _, _ = \
            statsmodels.stats.multitest.multipletests(pvals=self.p_value, alpha=q_crit, method='fdr_by')
        self.start_nucleotide = np.array(self.start_nucleotide)[significant]
        self.end_nucleotide = np.array(self.end_nucleotide)[significant]
        self.n_matched_bases = np.array(self.n_matched_bases)[significant]
        self.p_value = np.array(self.p_value)[significant]
        self.q_value = q_value[significant]

    def __len__(self):
        return len(self.start_nucleotide)


def compute_randomized_structures(transcript: Bio.Seq.Seq, n_randomizations: int, min_window: int,
                                  max_window: int, min_length: int,
                                  max_length: int, min_motif_length: int):
    raw_weights = np.array([transcript.count('A'), transcript.count('C'), transcript.count('G'), transcript.count('U')])
    weights = raw_weights / np.sum(raw_weights)

    def counting_wrapper(np_transcript, min_window, max_window, min_motif_length):
        motif_results = find_hairpins(np_transcript=np_transcript, min_window=min_window,
                                      max_window=max_window, min_motif_length=min_motif_length, parallel=False)
        if len(motif_results) == 0:
            return 0
        count = motif_results.count_matched_bases(start_motif=motif_results.order_start_position[0],
                                                  end_motif=motif_results.order_end_position[-1])
        return count

    def parallel_wrapper(seed):
        np.random.seed(seed)
        # Creates a random sequence according to the original distribution in the transcript
        random_sequence = np.random.choice([b'A', b'C', b'G', b'U'], size=max_length, p=weights)
        matched_base_count = [counting_wrapper(np_transcript=random_sequence[:length], min_window=min_window,
                                               max_window=max_window, min_motif_length=min_motif_length)
                              for length in range(min_length, max_length + 1)]
        return np.array(matched_base_count)

    seeds = np.random.randint(low=0, high=2 ** 32 - 1, size=n_randomizations)
    res = joblib.Parallel(n_jobs=joblib.cpu_count())(joblib.delayed(parallel_wrapper)(seeds[i])
                                                     for i in range(n_randomizations))
    # res = [parallel_wrapper(seeds[i]) for i in range(n_randomizations)]
    res_combined = np.vstack(res)
    return RandomizationResult(result_matrix=res_combined, min_length=min_length, max_length=max_length)


def predict_RNA_genes(sequence, n_randomizations: int, min_window: int, max_window: int, min_length: int,
                      max_length: int, q_crit: float,
                      min_motif_length) \
        -> IntervalResult:
    sequence = sequence.upper()
    transcript = sequence.transcribe()
    randomized_results = compute_randomized_structures(transcript=transcript, n_randomizations=n_randomizations,
                                                       min_window=min_window,
                                                       max_window=max_window,
                                                       min_length=min_length, max_length=max_length,
                                                       min_motif_length=min_motif_length)
    np_transcript = np.array(transcript, dtype=np.string_)
    actual_results = find_hairpins(np_transcript, min_window, max_window, min_motif_length)

    return randomized_results.compare(actual_results, q_crit=q_crit, len_sequence=len(sequence), min_window=min_window,
                                      min_motif_length=min_motif_length)


def write_results(results: IntervalResult, chrom: str, chromStart: int, res_prefix: str) -> None:
    # Creating *.bed file
    bed_name = res_prefix + '_track' + '.bed' if res_prefix != '' else 'track.bed'
    results.p_value
    with open(bed_name, mode='w') as f:
        header_line = '''track\tname=hairpins_results\tdescription="RNA genes predicted by hairpin"\tuseName=1\n'''
        f.write(header_line)
        for i in range(len(results)):
            this_chromStart = results.start_nucleotide[i] + chromStart
            this_chromEnd = results.end_nucleotide[i] + chromStart
            this_q_value = results.q_value[i]
            this_n_matched_bases = results.n_matched_bases[i]
            bed_line = '''{0}\t{1}\t{2}\t"q-value={3:.1e}"\n'''.format(chrom, this_chromStart, this_chromEnd,
                                                                       this_q_value)
            f.write(bed_line)
    rand_name = res_prefix + '_randomized_results' + '.txt' if res_prefix != '' else 'randomized_results.txt'
    with open(rand_name, mode='w') as f:
        header_line = "Length_of_interval\tmean_matched_bases\tstd_matched_bases\n"
        f.write(header_line)
        for i in range(len(results.randomized_result)):
            rand_line = "{0}\t{1}\t{2}\n".format(i + results.randomized_result.min_length,
                                                 results.randomized_result.mean_matched_bases[i],
                                                 results.randomized_result.std_error_matched_bases[i])
            f.write(rand_line)
    plot_name = res_prefix + '_plot' + '.pdf' if res_prefix != '' else 'plot.pdf'
    plt.figure()
    plt.errorbar(x=range(results.randomized_result.min_length, results.randomized_result.max_length + 1),
                 y=results.randomized_result.mean_matched_bases,
                 yerr=results.randomized_result.std_error_matched_bases, label='Randomized sequences')
    plt.plot(results.end_nucleotide - results.start_nucleotide,
             results.n_matched_bases, 'or', label='Real intervals reported')
    plt.xlabel('Length of interval')
    plt.ylabel('Number of matched bases in interval')
    plt.title('Comparison of the significant intervals found to the randomized results')
    plt.legend(loc='lower right')
    plt.savefig(plot_name, format='pdf')


@click.command()
@click.argument(u'infile')
@click.option('--n_randomizations', default=50, help="Number of randomizations used to assess the significance of the"
                                                     " results, corresponds to n_rand in the article",
              show_default=True)
@click.option('--min_window', default=4, help='Minimal window searching for matching nucleotides, corresponds to '
                                              'd_min in the article', show_default=True)
@click.option('--max_window', default=100, help='Maximal window searching for matching nucleotides, '
                                                'corresponds to d_max in the article', show_default=True)
@click.option('--min_length', default=20, help='Minimal length of reported RNA coding regions,'
                                               ' corresponds to r_min in the article', show_default=True)
@click.option('--max_length', default=100, help='Maximal length of reported RNA coding regions,'
                                                'corresponds to r_max in the article', show_default=True)
@click.option('--min_motif_length', default=4,
              help='The minimum number of consecutive base parings necessary to take such a motif into '
                   'consideration, corresponds to l_min in the article', show_default=True)
@click.option('--chrom', default='chr', help='String representation of chromosome name (such as chr10)',
              show_default=True)
@click.option('--chromStart', default=0, help='Start position of the sequence at the chromosome', show_default=True)
@click.option('--q_crit', default=0.05, help='Critical Benjamini–Yekutieli corrected false discovery '
                                             'rate for the results', show_default=True)
@click.option('--res_prefix', default='', help='Prefix for the output files', show_default=True)
def main(infile, min_window=4, n_randomizations=50, max_window=100, min_length=20, max_length=100, q_crit=0.05,
         min_motif_length=4,
         chrom='chr', chromstart=0, res_prefix=''):
    """
    :param INFILE: The filename (FASTA format) of the coding DNA strand being used to predict the RNA genes.
    """
    sequence = Bio.SeqIO.read(infile, 'fasta', alphabet=Bio.Alphabet.generic_dna).seq
    results = predict_RNA_genes(sequence, n_randomizations=n_randomizations, min_window=min_window,
                                max_window=max_window, min_length=min_length,
                                max_length=max_length, q_crit=q_crit, min_motif_length=min_motif_length)
    write_results(results=results, chrom=chrom, chromStart=chromstart, res_prefix=res_prefix)


if __name__ == '__main__':
    if sys.version_info.major == 3 and sys.version_info.minor >= 5:
        main()
    else:
        print('Please use Python 3.5 or later')
