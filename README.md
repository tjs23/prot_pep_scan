# prot_pep_scan
Scan protein sequences for permutations of amino acid subsequences


#proteome_consec_pep_scan.py

Finds sequences with the longest subsequences that match consecutive arrays of
sub-peptide sequences. You can run this script with something like:

`python3 proteome_consec_pep_scan.py /data/p3_ortho/seqs_metazoa/*.fasta.gz -p PEW PLP IRP GGP GPP -n 7 -o eukarya_hits.tsv`

There's some command-line too:

`$> python3 proteome_consec_pep_scan.py -h

usage: proteome_consec_pep_scan.py [-h] [-p PEP_SEQ [PEP_SEQ ...]]
                                   [-n PEP_COUNT] [-o OUT_TSV_FILE]
                                   FASTA_FILE [FASTA_FILE ...]


positional arguments:
  FASTA_FILE            One or more FASTA format sequence file paths
                        (separated by spaces). Wildcards accepted.

optional arguments:
  -h, --help            show this help message and exit
  -p PEP_SEQ [PEP_SEQ ...]
                        Peptide amino acid sub-sequences to search for. May
                        include "X" to match any amino acid. Space-separated,
                        without quotes. For example: RKL PGG STQ
  -n PEP_COUNT, --min-num-consec PEP_COUNT
                        Minimum number of consecutive sub-peptide sequences.
                        Default 3
  -g GAP_COUNT, --max-num-gaps GAP_COUNT
                        Maximum number of unspecified, internal sub-peptides;
                        gaps of sub-peptide length. Default 0
  -o OUT_TSV_FILE, --out-file OUT_TSV_FILE
                        Optional output file to write results as a tab-
                        separated table


For speed this uses the numba Python module (https://numba.pydata.org/).
