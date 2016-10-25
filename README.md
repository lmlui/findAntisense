# findAntisense
Find antisense RNA interactions

Authors: Lauren Lui, Andrew Uzilov, Andrea Corredor

Instructions for obtaining a platform and architecture specific c_antisense.so and be able to use findAntisense.py

1. Download the tar ball containing all source files
2. In the directory containing all source files type the command:
    python setup.py build_ext --inplace
3. Run findAntisense.py as described in findAntisense.py -h 

Description:

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

KNOWN BUG: Bad things if required sites exist outside of a duplex.
To fix (TODO):
Apply mask BEFORE generating duplex!  Required sites may exist OUTSIDE
of duplex (NOT even sub-duplex -- the PARENT duplex) in case of guide
overhang.  So, required sites may be in overhang!
