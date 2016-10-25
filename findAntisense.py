#!/usr/bin/env python2.7
"""
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
"""

# TODO:
# - Validate that mask string is correct.
# - Add way to load query seq from file (e.g. FASTA file).


import itertools
import operator
import optparse
import sys
import c_antisense as antisense
import feature
import os


# Base-pair conventions.
# TODO: These should be included from utils.h via SWIG somehow,
# otherwise they may go out-of-sync with this module.
bpCharNone = '.'
bpCharGu   = ':'
bpCharWc   = '|'

#DNA Alphabet set
DNA_alphabet = set('ACTG')

class Duplex():
    """Represents a sense/anti-sense interaction.
    """
    
    def __init__ (self, targetIndex, queryIndex, targetStart, queryStart,
                  duplexLen, annotBp, score):
        """Constructor.

        The object attributes exactly correspond to fields in the C layer
        "Duplex" struct.
        """
        self.targetIndex = targetIndex  # Index into target sequence list.
        self.queryIndex  = queryIndex   # Index into query sequence list.
        self.targetStart = targetStart  # Start of duplex in target RNA.
        self.queryStart  = queryStart   # Start of duplex in query RNA;
                                        # note this counts from 3' end of query!!
        self.duplexLen   = duplexLen    # Length of duplex.
        self.annotBp     = annotBp      # Base-pairing annotation string.
        self.score       = score        # Score of duplex.    


def loadFasta (fastaFileName, reqUniq=False):
    '''Read a FASTA file.

    Returns dict of name (string) -> sequence (string).
    '''
    print >> sys.stderr, 'Loading FASTA file "%s"...' % fastaFileName
    ret = dict()
    for line in open (fastaFileName, 'r'):
        if line[0] == '>':
            fields = line[1:].split()
            currentId = fields[0]
            if ret.has_key (currentId):
                if reqUniq:
                    print >> sys.stderr, \
                        'ID "%s" in FASTA input file "%s" not unique.' \
                        % (currentId, fastaFileName)
                    sys.exit (1)  # TODO: necessary?
                else:
                    print >> sys.stderr, \
                        'ID "%s" in FASTA input file "%s" not unique; overwriting previous sequence.' \
                        % (currentId, fastaFileName)
                    ret[currentId] = ''
            else:
                ret[currentId] = ''
        else:
            ret[currentId] += line.strip()
    return ret


def yieldFasta (fastaFileName):
    '''Read a FASTA file. Yields a tuple with (id,sequence).
    '''
    seq_list = [] #List tht holds sequence lines as they are read 
    first = True
    for line in open (fastaFileName, 'r'):
        if line.isspace() or line.startswith(';'): continue
        if line[0] == '>':
            if not first:
                yield (currentId,''.join(seq_list))
                seq_list = []
            fields = line[1:].split()
            currentId = fields[0]
            first = False
        else:
            seq_list.append(line.strip())

    #Handle last sequence
    yield (currentId,''.join(seq_list))
   


def score_alignment(align,wc,gu,o,e):
    """Scores alignment with affine model based on costs & gains passed."""
    score = 0
    in_gap = False
    for c in align:
        if c == '|':
            score += wc
            if in_gap: in_gap = False
        elif c == ':':
            score += gu
            if in_gap: in_gap = False
        elif c == '.' and in_gap: score -= e
        elif c == '.' and not in_gap:
            score -= o
            in_gap = True
    return score


# Create the command-line invocation parser.

def parse_arguments():
    """Create the command-line invocation parser."""
    parser = optparse.OptionParser (usage="""
        %prog [options] querySeq database.fasta database.bed
        """ + __doc__)

    parser.add_option ('--minLen', '-l', action='store', dest='minLen',type='int', metavar='LENGTH', default=10,help='minimum duplex length')

    parser.add_option ('--maxMis', '-m', action='store', dest='maxNumMismatches',type='int', metavar='LENGTH', default=1,help='maximum number of mismatches in duplex')

    parser.add_option ('--maxGU', '-g', action='store', dest='maxNumGUBasePairs',type='int', metavar='LENGTH', default=2,help='maximum number of G/U base pairs duplex')

    parser.add_option ('--mask', '-a', action='store', dest='mask',metavar='\'STRING\'', default=None,help='use this custom mask (put it in SINGLE quotes); default mask\
        is none, unless --cdMode is set, in which case a C/D sRNA-specific mask is created to constrain base pairing around the methylation site (as inferred by the N+5 rule)')

    parser.add_option ('--cdMode', '-s', action='store_true', dest='cdMode', default=False,help='turn on C/D sRNA-specific default mask creation and output options; if \
        this is set, query sequence is assumed to be the guide sequence that begins with the FIRST nuc of the D or D\ box (generally a "C")')

    parser.add_option ('--threshold', '-t', action='store', dest='threshold',type='int',default=10,help='Cutoff score used to filter duplexes. Default values is 10.')

    parser.add_option('--wc_gain','-w',help='Score gain for a single Watson-Crick pairing. Default is 2.',\
        default=2,type='float',dest='wc')

    parser.add_option('--gu_gain','-u',help='Score gain for a GU pairing. Default is 1.',default=1,\
        type='float',dest='gu')

    parser.add_option('--open','-o',help='Gap opening penalty. Default is 2.',default=2,type='float',dest='o')
    parser.add_option('--extend','-e',help='Gap extension penalty. Default is 1.',default=1,type='float',dest='e')
    parser.add_option('--to_stdout','-i',action='store_true',default=False,help='Output results to stdout.')
    return parser.parse_args()



def find_antisense(querySeqID,querySeq,fastaFileName,targetFileName,minLen,maxNumMismatches,maxNumGUBasePairs,\
		mask,cdMode,threshold,wc,gu,o,e,to_stdout):
    """Finds antisense targets. Arguments correspond directly to the command line arguments."""
    # Make query into RNA.
    querySeq = querySeq.upper().replace ('T', 'U')  # TODO: should be a lib func

    # If user provided a custom mask, make sure it makes sense.
    if mask is not None:
        if (len (mask) != len (querySeq)):
            raise Exception ("Mask and query sequence must be of same length.")

    # Echo search criteria, for logging.
    print >> sys.stderr, \
        'Searching through file "%s" for anti-sense match to seq %s' \
        % (fastaFileName, querySeq)
    print >> sys.stderr, \
        'Minimum duplex length is set to %d' % minLen
    print >> sys.stderr, \
        'Maximum number of mismatches is set to %d' % maxNumMismatches
    print >> sys.stderr, \
        'Maximum number of GU base pairs is set to %d' % maxNumGUBasePairs
    if cdMode:
        print >> sys.stderr, "Running in C/D sRNA mode."
    if mask is not None:
        print >> sys.stderr, \
            'User specified custom mask for query seq:\n%s\n%s' \
            % (querySeq, mask)

    if mask is not None:
        # Use user-specified mask.
        print >> sys.stderr, "Using user-specified mask:\n%s\n%s" \
            % (querySeq, mask)
    elif cdMode:
        # Make a mask for this seq that disallows mismatches and G/U
        # base pairs at methylation site as inferred by N+5 rule.
        #
        # IMPORTANT: We always assume that last char of 'seq' string is the 1st char
        # of the D or D' box!
        assert len (querySeq) >= 6, \
            "Can't make mask for a guide (query) of %d nt because can't apply N+5 rule for something < 6 nt." \
            % len (querySeq)

        mask = ['.' for i in xrange (len (querySeq))]
        if len (mask) >= 7:
            mask[-7] = '!'  # next (5') to meth site
        mask[-6] = '!'      # 2'-O-me site
        mask[-5] = '!'      # next (3') to meth site
        mask = ''.join (mask)
        print >> sys.stderr, "Made mask using N+5 rule:\n%s\n%s" \
            % (querySeq, mask)
    else:
        # Make a default, completely permissive mask (effectively no mask).
        mask = ''.join (['.' for i in xrange (len (querySeq))])
        print >> sys.stderr, "Made default no-op mask:\n%s\n%s" \
            % (querySeq, mask)

    # Set us up the config.
    # Input to antisene.setConfig() MUST be two lists of strings.
    # The 1st list contains option names, the 2nd option values.
    # Lists are in 1-to-1 correspondence.
    # Option values MUST be strings, even if they are numbers.
    #
    # NB: We reverse the mask (make it 3'->5') because it corresponds to nucs
    # in the query, and the query gets reversed later.

    ret = antisense.setConfig (
        ['mode', 'minLen', 'maxGuBp', 'maxMismatches', 'mask'],
        [str (x) for x in (
            1, minLen, maxNumGUBasePairs, maxNumMismatches,
            mask[::-1]
        )])
    print >> sys.stderr, 'Set up the config, returned error:'
    print >> sys.stderr, ret
    assert ret is None  # don't proceed if error was returned


    # Load targets.
    print >> sys.stderr, 'Loading targets from "%s"...' % targetFileName
    targets = feature.bedToLocusColl (
        targetFileName, useId=True)
    for id, seq in loadFasta (fastaFileName).iteritems():
        targetSeq = seq.upper().replace ('T', 'U')
        targets.getFeatureByID (id).seq = targetSeq


    # Find sense/anti-sense duplexes between qery and targets that are at least
    # as good as criteria we specified in antisense.setConfig().
    #
    # NB: Before calling the C function antisense.duplexes(), we reverse the query
    # (make it 3'->5') because the C func requires that to slide it against
    # 5'->3' target seqs.

    print >> sys.stderr, 'Searching for targets to query "%s"...' % querySeq

    duplexes = []
    for d in antisense.findDuplexes ([t.seq for t in targets.featuresOrdered],
                                     [querySeq[::-1]]):
        # Convert tuple "d" (created by SWIG from the C linked list "DuplexList")
        # to the Python "Duplex" object.
        #
        # TODO: Ideally, this conversion should happen in the SWIG layer.
        # Unfortunately, Andrew doesn't know how to do that (yet), so instead he
        # makes simple Python tuples in the SWIG layer (which is easier) and uses
        # those to ferry duplex data from the C layer to the Python layer.
        
        # Index into query list should always be 0
        # (only 1 query, list contains 1 item).
        assert d[1] == 0

        # Make Duplex object.
        # See the Duplex class for explanation of how these work.
        duplexes.append (Duplex (targetIndex = d[0],
                                 queryIndex  = d[1],
                                 targetStart = d[2],
                                 queryStart  = d[3],
                                 duplexLen   = d[4],
                                 annotBp     = d[5],
                                 score       = d[6]))

    print >> sys.stderr, 'Computed %d duplexes in C layer...' % len (duplexes)


    # Rescore and filter duplexes.
    #
    # TODO: This is for prototyping.  Once a good scoring/filtering scheme is
    # figured out, it should be ported to C layer and happen there.

    filteredDuplexes = []
    for duplex in duplexes:
        duplex.score = score_alignment(duplex.annotBp,wc,gu,o,e)
        # TODO: FILTERING HERE.
        if duplex.score > threshold: filteredDuplexes.append(duplex)

    print >> sys.stderr, \
          'Filtered %d duplexes down to %d duplexes in Python layer...' \
          % (len (duplexes), len (filteredDuplexes))


    # Sort filtered duplexes by their scores.
    # Technically speaking, we sort the list in place, comparing on the values of
    # the 'score' attribute of the objects in the list.
    filteredDuplexes.sort (key=operator.attrgetter ('score'),reverse=True)


    # Print filtered, ordered duplexes.
    if to_stdout: f = sys.stdout
    else:
        f = open(''.join([querySeqID,'.results.txt']),'w')

    for duplex in filteredDuplexes:
        # This is the ENTIRE target RNA, not just the part that is in the duplex.
        target = targets.featuresOrdered[duplex.targetIndex]
        
        # Print the duplex.
        #
        # The key thing to remember about layout is that the duplex display width
        # is determined solely by query length.
        # We are going to display the entire query sequence no matter what, and we
        # will fit as much target sequence into that space as we can.
        #
        # All output coords (except BED line) are in 1-based inclusive coords.
        # Internally, all coords are 0-based half-open (in-between).

        if cdMode:
            # Line showing methylation position.
            methPosLine  = '     *'
            methPosLine += ''.join (
                itertools.repeat (' ', len (querySeq) - len (methPosLine)))
            methPosLine += "  2'-O-me site"
            print >> f, methPosLine

        # Compute length of query sequence outside of duplex, i.e. "extra" sequence
        # to be displayed on left (3') and right (5') of query.
        extraL = duplex.queryStart
        extraR = len (querySeq) - (duplex.queryStart + duplex.duplexLen)
        assert extraL + duplex.duplexLen + extraR == len (querySeq)

        # Target sequence to display.
        # Display as much target sequence as query sequence, being careful not go
        # go out of bounds.
        targetL = max (0, duplex.targetStart - extraL)
        targetR = min (
            len (target), duplex.targetStart + duplex.duplexLen + extraR)
        targetSeq = target.seq[targetL:targetR]
      
        # If query overhangs end(s) of target sequence, we need to pad target
        # so it aligns to query in the output. 
        targetPadL = max (0, -(duplex.targetStart - extraL))
        targetPadR = max (
            0, duplex.targetStart + duplex.duplexLen + extraR - len (target))
     
        targetLine = "%s%s%s  5'->3' target \"%s\" (duplex: [%d,%d])" % (
            ''.join (itertools.repeat (' ', targetPadL)),
            targetSeq,
            ''.join (itertools.repeat (' ', targetPadR)),
            target.id,
            duplex.targetStart + 1,
            duplex.targetStart + duplex.duplexLen)
        print >> f, targetLine
        
        bpLine = '%s%s%s  duplexLen=%d WC=%d GU=%d mismatches=%d' % (
            ''.join (itertools.repeat (' ', extraL)),  # left padding
            duplex.annotBp,
            ''.join (itertools.repeat (' ', extraR)),  # right padding
            duplex.duplexLen,
            duplex.annotBp.count (bpCharWc),
            duplex.annotBp.count (bpCharGu),
            duplex.annotBp.count (bpCharNone))
        print >> f, bpLine
        
        queryLine = "%s  3'->5' %s (duplex: [%d,%d])" % (
            querySeq[::-1],
            ('guide' if cdMode else 'query'),
            len (querySeq) - (duplex.queryStart + duplex.duplexLen) + 1,
            len (querySeq) - duplex.queryStart)
        print >> f, queryLine

        maskLine = mask[::-1] + '  mask'
        print >> f, maskLine

        # Counts distance from D or D' box, but only up to 9 nt away
        # (since we start running out of single-digit numbers at that point).
        countFromBoxLine  = ''.join ([
            str (i) for i in xrange (0, min (10, len (querySeq)))])
        countFromBoxLine += ''.join (
            itertools.repeat (' ', len (querySeq) - len (countFromBoxLine)))
        if cdMode:
            countFromBoxLine += "  distance from 1st nuc of D or D' box"
        else:
            countFromBoxLine += "  distance from 3'-most nuc of query"
        print >> f, countFromBoxLine
        
        # BED line for the target region (duplex part only) in genomic coords.
        targetDuplexOnly = feature.LocatedFeature (
            target,
            duplex.targetStart,                     # start rel to target
            duplex.targetStart + duplex.duplexLen,  # end rel to target
            coordSystem = 2,
            strand = '+',  # subfeat is on plus strand of its parent
            name = target.id,
            score = duplex.score)
        targetDuplexOnly.mapToParentRef()  # convert local to genomic coords
        print >> f, "BED of target (duplex part only):"
        print >> f, targetDuplexOnly.toBed (numFields=6, nameField='name')

        # Stanza separator for good measure.
        print >> f, "//"

    if to_stdout:
        #print separator for clarity
        print >> f, '######################################################################################'
    else:
        f.close()

def main(args):
    """Main function."""
    # Parse command-line invocation.
    opts, args = parse_arguments()
    if len (args) != 3:
        raise Exception(
    'Incorrect number of required positional arguments; need 3 arguments, you provided %d and they were:\n    %s'
    % (len (args), ', '.join (args)))
    querySeq       = args[0]
    fastaFileName  = args[1]
    targetFileName = args[2]

    if os.path.isfile(querySeq): #If argument is a fasta file
        assert opts.mask is None, "Can't pass a mask if guide sequences being read from a file."
        with open(querySeq,'r') as query_file:
            print >> sys.stderr, "Iterating over guide sequences in file %s" % querySeq
            for seq in yieldFasta(querySeq):
                find_antisense(seq[0],seq[1],fastaFileName,targetFileName,opts.minLen,opts.maxNumMismatches,\
                                opts.maxNumGUBasePairs,opts.mask,opts.cdMode,opts.threshold,opts.wc,opts.gu,opts.o,opts.e,opts.to_stdout)
		#separator for error messges
		print >> sys.stderr, '************************************************'


    elif all(c.upper() in DNA_alphabet for c in querySeq): #Argument is a single sequence
	#No guide sequence id
	querySeqID = querySeq
        find_antisense(querySeqID,querySeq,fastaFileName,targetFileName,opts.minLen,opts.maxNumMismatches,\
        opts.maxNumGUBasePairs,opts.mask,opts.cdMode,opts.threshold,opts.wc,opts.gu,opts.o,opts.e,opts.to_stdout)

    else: 
        #String passed is neither a file nor a valid DNA sequence.  
        #Proceed as though it were a sequence, but inform user of anomaly
        print >> sys.stderr, "Guide sequence passed %s does not look like DNA, but no file found with that name. Proceeding anyway." % querySeq
	#No guide sequence id
	querySeqID = querySeq
        find_antisense(querySeqID,querySeq,fastaFileName,targetFileName,opts.minLen,opts.maxNumMismatches,\
        opts.maxNumGUBasePairs,opts.mask,opts.cdMode,opts.threshold,opts.wc,opts.gu,opts.o,opts.e,opts.to_stdout)


if __name__ == "__main__" :

    try:
        sys.exit(main(sys.argv))
    except EnvironmentError as (errno,strerr):
        sys.stderr.write("ERROR: " + strerr + "\n")
        sys.exit(errno)


