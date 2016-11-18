# findAntisense
Find antisense RNA interactions!

Authors: Lauren Lui, Andrew Uzilov, Andrea Corredor

**Note that this will only work with Python 2.7**

**Installation**

Instructions for obtaining a platform and architecture specific c_antisense.so and be able to use findAntisense.py

1. Download the zip file containing all source files
2. In the directory containing all source files type the command:
    python setup.py build_ext --inplace
3. Run findAntisense.py as described in findAntisense.py -h 

**Examples and use:**

See the wiki: https://github.com/lmlui/findAntisense/wiki

**Description:**

Given a single sequence as a query, search a sequence database for
sense/anti-sense interactions (duplexes) between query and target(s).
Currently, the lengths of sense and anti-sense sequences in predicted duplexes
are the same (no bulges), but this may change in the future.
All input sequences are assumed to be written 5' to 3' left to right.  We do
not search reverse complements of any given sequence, so make sure to give the
correct strand.

We display sequences as RNA, but in theory the searching should work with DNA,
as long as the specified config (thresholds, match criteria, scoring function,
etc.) makes sense for whatever you are searching.

The BED file should contain coordinates for sequences in the FASTA file so we
can compute and output coordinates for matches to targets.  The names in the
FASTA headers must match the names in the BED file so we can associate
sequences with their coordinates.

Masks can be specified to allow/disallow base-pairing at specific query
positions.  The mask string must be single-quoted to prevent shell expansion
and must be of same length as the query seq.  It must be made up of only the
following chars that denote constraints at each position:

    ! = only W/C base pairing allowed
    
    ? = any canonical base pairing (W/C or G/U) allowed
    
    . = no constraints (anything allowed)
    
Duplexes that pass the match criteria are written to STDOUT.  Search criteria
are written to STDERR for logging.

Output is by default written to file(s) <seqID>.output.txt when the guide 
sequences are read from a FASTA file. When a single guide sequence is 
passed--and no sequence id is known--the output file will be <guide sequence>.output.txt
Results can also be outputted to stdout by passing the corresponding flag. 

Options:

    -h,        --help             show this help message and exit 
    -l LENGTH, --minLen=LENGTH    minimum duplex length
    -m LENGTH, --maxMis=LENGTH    maximum number of mismatches in duplex
    -g LENGTH, --maxGU=LENGTH     maximum number of G/U base pairs duplex
    -a 'STRING', --mask='STRING'  use this custom mask (put it in SINGLE quotes);
                                  default mask is none, unless --cdMode is set,
                                  in which case a C/D sRNA-specific mask is created to
                                  constrain base pairing around the methylation site (as
                                  inferred by the N+5 rule)
                        
    -s, --cdMode                  turn on C/D sRNA-specific default mask creation and
                                  output options; if this is set, query sequence
                                  is assumed to be the guide sequence that begins with
                                  the FIRST nuc of the D or D\ box (generally a "C")
                        
    -t THRESHOLD, --threshold=THRESHOLD
                                  Cutoff score used to filter duplexes. Default values
                                  is 10.
                        
    -w WC, --wc_gain=WC           Score gain for a single Watson-Crick pairing. Default
                                  is 2.
                        
    -u GU, --gu_gain=GU           Score gain for a GU pairing. Default is 1.
  
    -o O, --open=O                Gap opening penalty. Default is 2.
    -e E, --extend=E              Gap extension penalty. Default is 1.
    -i, --to_stdout               Output results to stdout.


KNOWN BUG: Bad things if required sites exist outside of a duplex.
To fix (TODO):
Apply mask BEFORE generating duplex!  Required sites may exist OUTSIDE
of duplex (NOT even sub-duplex -- the PARENT duplex) in case of guide
overhang.  So, required sites may be in overhang!
