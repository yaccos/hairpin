# Purpose
Hairpin is a program which predicts RNA genes in a genome _ab-intio_.
It is based on the assumption that regions in the genome having more possible
pairings than expected by chance, will possible code for structured RNA genes.
The resulting `*.bed`-file can be uploaded as a track to a genome browser based
upon the same genome assembly.

# Requirements

This is intended to work for Python 3.5.0 or later. 
The required packages with their versions are listed in `requirements.txt`,
run `pip3 install -r requirements.txt` in order to install them.

# Usage
The program is run from the command line using the syntax (example):
```bash
python3 hairpin.py --n_randomizations 50 --min_window 4 --max_window 20 \
--min_length 20 --max_length 100 --min_motif_length 4 --chrom 'chr' \
--chromStart 1000 --q_crit 0.05 --res_prefix 'my_results' e_coli_genome.fasta
```
The arguments are described upon calling:

```bash
python3 hairpin.py --help
```

Also, you *can* use the functions and classes provided by the software as part of a custom script as done in the
examples. Note though that the internal API is poorly documented, so issues and pull requests regarding its documentation
are welcome.

## More about the arguments
Except the input file, all options have default values. 
The input must be in `fasta`-format, 
where the sequence shows the coding DNA strand on a chromosome. An input file with multiple records is allowed,
but only the first record will be analyzed. Nucleotide ambiguity codes are accepted, but ignored in the analyses.
In order to provide specifications to the genome browser,
the name of the chromosome or scaffold must be specified much as `--chrom chr10`. The default argument however,
is OK for procaryotes only having a single genome. In case the record is not the entire chromosome, 
consider the following cases:

* The record does cover the start of the chromosome: Specify the start position (zero-indexed) as
 (for instance) `--chromStart 999`
* The record does not cover the end of the chromosome: No action required

## Note on performance

This program uses parallelization, meaning it might take up huge amounts of the computer's working capacity. For large 
input files or many rounds of randomizations, running on a computer with many cores is recommended. 

# Output
The program outputs three files:

* `(res_prefix)track.bed`: The genomic track information, can be uploaded to a genomic track browser
* `(res_prefix)randomized_result.txt`: Tab delimited text file specifying the results when constructing
randomized sequences
* `(res_prefix)plot.pdf`: Plot showing the number of matched bases for randomizations and
 reported significant intervals. The blue line correspond to the mean number of matched bases for the
  randomizations, while the error bar denotes its standard deviations. Each red dot corresponds to a
  significant interval found.
  
# Disclaimer
Do not expect the program to yield any meaningful biological results. As the referenced paper states, 
the biological signal is typically too small to be separated from random noise. You are welcome to prove the author
wrong though.


# References
Hairpin: a tool for predicting structural non-coding RNA genes _ab initio_
Jakob Peder Pettersen
bioRxiv 640862; doi: https://doi.org/10.1101/640862
