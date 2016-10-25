#!/bin/env python
"""
Input/output in various bioinformatics data formats.

For reading data, this class should provide generators that read the data and
yield standard Python data types.  So, DO NOT yield objects of custom classes
or write factory functions that return custom collections.
The philosophy is: provide generators that yield simple, standard data
structures and leave it up to the caller to figure out how to (or even
whether to) use them to instantiate and/or populate whatever objects they want.
The readers here are to simplify the READING itself, not to assemble data into
some specific data model (which is dangerous anyway, because one should not
guess at what the caller wants to do with the data).
"""

# TODO:
# - Option for writing to gzipped file directly (wiggles can get large, but
#   compress real easy).
# - Debug countdown goes backwards when doing minus strand.
# - Many asserts should be more gentle (e.g. raise exceptions useful to user
#   instead).  Same with sys.exit() calls.
# - Standardize naming (fix stuff like: iterFasta vs. loadWiggleAsIter).


import gzip
import itertools
import re
import sys

import exceptions
import feature
from misc import printe


logPrefix = '[binfio]'


###############################################################################
# Helpers.
###############################################################################


def openForReading (fileName):
    """Open a file for reading, trying to figure out whether it is a gzip, bz2,
    or plain file, and return its file object, so the caller does not have to
    bother figuring out whether the file is compressed or uncompressed.

    The file objects for gzip and bz2 files support most, but not all, of the
    methods of plain files (created by the built-in 'open' function).  See the
    docs for 'gzip' and 'bz2' modules for details on that.  But for the most
    part, they can be treated like regular file objects.

    Universal newline support works for all files except gzip.
    """
    try:
        f = open (fileName, 'r')
        # If file does not exist or is not readable, an exception will be
        # raised at this point.
    except IOError, e:
        raise e
    f.close()

    # If we made it this far, file exists and is readable.
    # Try to figure out the file type, open the file, and return the file
    # object.
    
    try:
        # gzip
        f = gzip.open (fileName, 'rb')  # TODO: why goes 'U' flag not work?
        f.read (1)  # try to read 1 byte; will fail for non-gzipped files
        f.close()
        return gzip.open (fileName, 'rb')
    except IOError:
        pass  # not a gzip file

    # bz2
    try:
        # Do the import here and catch if it fails, because the 'bz2' module
        # may not be installed.
        import bz2
        try:
            f = bz2.BZ2File (fileName, 'rbU')
            f.read (1)  # try to read 1 byte; will fail for non-bz2 files
            f.close()
            return bz2.BZ2File (fileName, 'rbU')
        except IOError:
            pass  # not a bz2 file
    except ImportError, e:
        printe ('WARNING: "bz2" module cannot be imported (maybe not installed?)  Will not be able to open bz2 files.  Exception text:\n' + str (e))

    # If we made it this far, file must be a plain file.
    return open (fileName, 'rU')


###############################################################################
# BED.
# File format definition here:
# http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED
#
# TODO: This parser needs to be tested by a bunch of malformed BED files!!
# Not every single failure case has been tested, just that well-formed BED
# files tend to pass.
###############################################################################


class BedParseError (Exception):
    """Exception that indicates a malformed BED file.
    """
    def __init__ (self, lineNum, fileName, line, msg=None):
        self.lineNum  = lineNum
        self.fileName = fileName
        self.line     = line
        self.msg      = msg

    def __str__ (self):
        if self.fileName is None:
            inputName = 'BED input from STDIN'
        else:
            inputName = 'BED file "%s"' % self.fileName
            
        if self.msg is None:
            return 'Parse error occured on line %d of %s.  Offending line:\n%s' \
                   % (self.lineNum, inputName, self.line)
        else:
            return 'Parse error occured on line %d of %s.  Offending line:\n%s\n%s' \
                   % (self.lineNum, inputName, self.line, self.msg)


# TODO: Add a 'disableParanoia' option for more efficient loading.  Will
# require quite a bit of refactoring.
def iterBed (fileName=None, numCols=6, trimWhitespace=True, ignoreStrand=False,
             verbose=0):
    """Parse BED data, with many safety checks.

    This is a generator that, for each BED line, yields a dictionary mapping
    BED field names (C{chrom, chromStart, chromEnd, name, score, strand,
    thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts}) to their
    values.

    If C{fileName} is provided, input is read from that file.  It can either be
    a string path to a file, in which case we try to open it, or an iterable
    yielding string BED lines (this allows the caller to e.g. open a file and
    just pass the open file object as C{fileName}).
    If a path is given, input can be be uncompressed or compressed (C{gzip} or
    C{bz2}) files. 
    The file type is figured out based on file contents at that path, not name.
    If C{fileName} is not provided, or it is C{None}, input is read from
    C{STDIN}, in which case it must be uncompressed.
    
    If C{numCols} is provided, only the first C{numCols} fields are loaded and
    validated; if the input data contains more than C{numCols} fields on some
    line, they are silently ignored.
    Also, each line must contain at least C{numCols} fields.
    If C{numCols} is C{None}, all fields on a line are loaded, and BED input can
    actually contain lines with variable numbers of fields; the actual number
    of fields on each line is saved into its dictionary under the key
    C{numFields}.
    Note that C{numFields} is set even if a C{numCols} number was provided, for
    consistency.
    
    All BED field mappings get a value, but ignored values are set to C{None}.
    This way, all keys corresponding to BED fields are defined, so the caller
    does not have to check for C{KeyError} exceptions.
    For example, if C{numCols} is 6, then keys C{thickStart} and beyond will map
    to C{None} values.
    
    Fields are assumed to be tab-delimited (NOT whitespace delimited).  As a
    check against various forms of manual editing insanity, we by default also
    trim the leading and trailing whitespace from each field, although that can
    be disabled by setting C{trimWhitespace} to C{False}, should you decide that
    e.g. identifiers containing leading/trailing whitespace are a good idea.

    Much type conversion and validation is done on BED field values.
    The following fields become ints: C{chromStart, chromEnd, thickStart,
    thickEnd, blockCount}.
    The BED score column is treated as a float and the range is not validated,
    (despite the specification stating it is an int in range C{[0,1000]})
    because people tend to squeeze all kinds of score types in there.
    Comma-separated fields (C{itemRgb, blockSizes, blockStarts}) are split by
    comma, converted to ints, validated, and each value becomes a list.

    Safety checks on the strand column are disabled if you set
    C{ignoreStrand} to True, but the string value in the column is still
    returned.  This allows you to have strandless features with intron/exon
    notation, or use a different strand convention.
    """
    if numCols is None:
        callerSetNumCols = False
    else:
        callerSetNumCols = True
        if numCols < 3:
            raise exceptions.InvalidArgumentError, \
                  "Cannot have less than 3 columns in a BED file (must have numCols >= 3, you said numCols=%d)." \
                  % numCols

    if fileName is None:
        # Read from STDIN.
        inFile   = sys.stdin
        fileName = 'STDIN'
    else:
        if isinstance (fileName, str):
            # Read from a non-STDIN file.
            inFile = openForReading (fileName)
        else:
            # Assume caller already gave us data as an iterable ready-to-go.
            inFile   = fileName
            fileName = '(no name)'

    lineNum = 0
    for line in inFile:
        lineNum += 1

        if verbose == 1:
            if lineNum % 1000000 == 0:
                print >> sys.stderr, \
                      logPrefix, 'Loading line %d of BED input from "%s"...' \
                      % (lineNum, fileName)
        elif verbose == 2:
            if lineNum % 1000 == 0:
                print >> sys.stderr, \
                      logPrefix, 'Loading line %d of BED input from "%s"...' \
                      % (lineNum, fileName)
        
        if line.startswith ('track') or line.startswith ('browser'):
            continue  # skip header lines
        if line.strip() == '':
            continue  # skip blank lines
        
        fields = line.strip().split ('\t')
        if callerSetNumCols:
            # Caller specified a number of columns to expect, so need to check
            # that input has at least that many.
            if len (fields) < numCols:
                raise BedParseError (
                    lineNum, fileName, line,
                    "Wrong number of fields; should be at least %d fields, found %d."
                    % (numCols, len (fields)))
        else:
            # Caller did not specify how many columns to expect.
            if len (fields) < 3:
                raise BedParseError (
                    lineNum, fileName, line,
                    "BED line contains only %d columns.  Must have a minimum of 3 columns." \
                    % (len (fields)))
            else:
                # Save the actual number of columns on this line, since the rest
                # of this function relies on it so much.
                numCols = len (fields)

        # Clean up string fields.
        if trimWhitespace:
            for f in (0, 3, 5):
                if f < numCols:
                    fields[f] = fields[f].strip()
                
        # Validate int fields.
        for f in (1, 2, 6, 7, 9):
             if f < numCols:
                 try:
                     fields[f] = int (fields[f])
                 except (TypeError, ValueError):
                     raise BedParseError (
                         lineNum, fileName, line,
                         'Expected integer in field %d, but found this: "%s".'
                         % (f+1, fields[f]))

        # Other obsessive validation.

        if fields[1] < 0:
            raise BedParseError (
                lineNum, fileName, line, "Field 2 (chromStart) must be >= 0.")

        if fields[2] < 0:
            raise BedParseError (
                lineNum, fileName, line, "Field 3 (chromEnd) must be >= 0")

        if fields[1] > fields[2]:
            raise BedParseError (
                lineNum, fileName, line,
                "Field 2 (chromStart) must be <= field 3 (chromEnd).")

        if numCols >= 5:
            try:
                fields[4] = float (fields[4])
            except (TypeError, ValueError):
                raise BedParseError (
                    lineNum, fileName, line,
                    'Expected float in field 5 (score), but found this: "%s".'
                    % (fields[4]))

        if numCols >= 6 and not ignoreStrand:
            if fields[5] != '+' and fields[5] != '-':
                raise BedParseError (
                    lineNum, fileName, line,
                    'Invalid field 6 (strand); must be either "+" or "-", but it is: "%s".' \
                    (fields[5]))

        if numCols >= 8:
            if fields[6] < fields[1] or fields[7] > fields[2]:
                raise BedParseError (
                    lineNum, fileName, line,
                    "Coords in columns 7 and 8 (thickStart, thickEnd = %d, %d) are out of feature bounds (chromStart, chromEnd = %d, %d)."
                    % (fields[6], fields[7], fields[1], fields[2]))

        if numCols >= 9:
            rgb = fields[8].strip().split (',')
            if len (rgb) == 1:
                # Apparently you can have a single RGB value sometimes.  I guess
                # it means "no RGB".  Go figure, the spec doesn't explicitly
                # mention it.
                # If this happens, hack it to be a black triplet.
                rgb = [0, 0, 0]
            else:
                # Parse and validate the RGB triplet.
                if len (rgb) != 3:
                    raise BedParseError (
                        lineNum, fileName, line, "Invalid field 9 (RGB).")
                try:
                    rgb = [int (x) for x in rgb]
                except (TypeError, ValueError):
                    raise BedParseError (
                        lineNum, fileName, line,
                        "Field 9 (RGB) does not contain all integers.")
                if min (rgb) < 0 or max (rgb) > 255:
                    raise BedParseError (
                        lineNum, fileName, line,
                        "Field 9 (RGB) contains value(s) out of [0,255] range.")
            fields[8] = rgb

        if numCols >= 10:
            blockCount  = fields[9]
            if blockCount < 1:
                raise BedParseError (
                    lineNum, fileName, line,
                    "Field 10 (blockCount) should be >= 1, but it is %d.  (For features with no introns, it should be 1 because the whole feature is an exon.)"
                    % (blockCount))
            fields[9] = blockCount

        if numCols >= 12:
            blockSizes  = fields[10].strip().split(',')
            blockStarts = fields[11].strip().split(',')

            # The Table Browser seems to leave a trailing empty list item
            # (e.g. '1414,2194,' for 2 blocks), so we must remove it if it
            # exists.
            if blockSizes[-1] == '':
                blockSizes.pop()
            if blockStarts[-1] == '':
                blockStarts.pop()
                
            if len (blockSizes) != blockCount or len (blockStarts) != blockCount:
                raise BedParseError (
                    lineNum, fileName, line,
                    "Number of blocks in fields 11 and/or 12 does not match block count in field 10.")
            try:
                blockSizes  = [int (x) for x in blockSizes]
                blockStarts = [int (x) for x in blockStarts]
            except (TypeError, ValueError):
                raise BedParseError (
                    lineNum, fileName, line,
                    "Fields 11 and 12 are not all integers.")

            # TODO: Validate that all blocks are within bounds.
            # TODO: Validate that blocks don't overlap.
            # TODO: Validate that all blocks are in order; if not, put them in
            # order for sanity.

            fields[10] = blockSizes
            fields[11] = blockStarts

        # Construct and return BED fields as a dictionary.

        fieldNames = ['chrom', 'chromStart', 'chromEnd', 'name', 'score',
                      'strand', 'thickStart', 'thickEnd', 'itemRgb',
                      'blockCount', 'blockSizes', 'blockStarts']
        fieldsDict = {}

        # Save values loaded from BED data.
        for i in xrange (numCols):
            fieldsDict[fieldNames[i]] = fields[i]

        # Put None-type placeholders for remaining values, so that caller does
        # not have to deal with KeyError exceptions.
        for j in xrange (numCols, len (fieldNames)):
            fieldsDict[fieldNames[j]] = None

        fieldsDict['numFields'] = numCols

        yield fieldsDict
        


def featToBed (feat, numFields, nameField=None, itemRgb=0):
    """Convert feature to a tab-delimited BED format string with C{numFields}
    fields (columns), without a newline terminator, which is then returned.

    The type of C{feat} is never checked; instead, we use duck typing and assume
    the caller passed in the appropriate object that contains attributes we need
    for the needed BED fields.  The appropriate names are those used by my usual
    featured classes (L{feature.LocatedFeature}, L{feature.Spliced}).

    Scores are converted to integers because the C{bedToBigBed} tool choked on
    floats at the time of this writing.
    If there is no score, the integer 1000 (max BED score) is output as a
    dummy placeholder.

    If C{itemRgb} is passed in, it is used for the 9th field.
    This input is NOT validated; the user must provide a correct value.
    Otherwise, it defaults to a C{0} placeholder.

    @todo: Controlling the following fields need some work (currently just set
        to some sensible default): thickStart, thickEnd.

    @param nameField: For BED lines with 4 fields or greater, this argument
        specifies how to get the BED name for the 4th field.
        If C{nameField} is a string set to either C{name} or C{id}, we use that
        feature property for the BED name.
        Otherwise, we use C{nameField} as a callable function that takes as input
        C{feat} and returns a string name.
        If C{nameField} is None, we will try the C{id} property, then C{name}
        property, then use a single dot as a placeholder if both are None.
    """
    if not (3 <= numFields <= 12) \
            or numFields == 7 or numFields == 10 or numFields == 11:
        raise exceptions.InvalidArgumentError, \
              "Number of BED fields (%d) is invalid; valid values are: 3, 4, 5, 6, 8, 9, 12." \
              % (numFields)
    
    retFields = [feat.landmarkName, str (feat.L), str (feat.R)]
    
    if numFields >= 4:
        # Need a name for the 4th column.
        if nameField is not None:
            # User told us how to get the name.
            if nameField == 'id' or nameField == 'name':
                name = getattr (feat, nameField)
            else:
                name = nameField (feat)
        else:
            # Need to figure out how to get name.
            # Try ID if one exists, then name if exists, then a
            # placeholder/dummy name.
            if feat.id is not None:
                name = feat.id
            elif feat.name is not None:
                name = feat.name
            else:
                # Placeholder.
                name = '.'
        retFields.append (name)

    if numFields >= 5:
        # Need a score for the 5th column.
        if feat.score is not None:
            # TODO: Ideally should do this assert, but that ties our
            # hands, so screw it for now.
            #assert feat.score >= 1000 and feat.score <= 0
            # BED scores require ints
            # (bedToBigBed chokes on non-ints, probably other tools do, too).
            score = str (int (feat.score))
        else:
            # Dummy score.
            score = '1000'
        retFields.append (score)

    if numFields >= 6:
        # Need a strand for the 6th column.
        if feat.strand is not None:
            strand = feat.strand
        else:
            # Placeholder.
            strand = '.'                        
        retFields.append (strand)

    if numFields >= 8:
        # Default thickStart/thickEnd to start/end of feature (TODO).
        retFields.append (str (feat.L))
        retFields.append (str (feat.R))

    if numFields >= 9:
        retFields.append (str (itemRgb))

    if numFields >= 12:
        retFields.append (str (len (feat.exons)))  # blockCount
        if feat.strand == '-':
            # Spliced feats store exons in order from 5' to 3' end of the
            # feature, but BED files store from lower to higher chromosome
            # coordinate regardless of strand.
            # So, need to change exon blocks to be relative to lower coord
            # (3' end of feature).
            exons = feature.reverseExons (feat.exons, len (feat))
        else:
            exons = feat.exons
        # blockSizes
        retFields.append (','.join ([str (e[1] - e[0]) for e in exons]))
        # blockStarts
        retFields.append (','.join ([str (e[0]) for e in exons]))

    return '\t'.join (retFields)


def collToBedFile (coll, fileName=None, compress=None, numFields=6,
                   nameField=None):
    """Write a feature collection to a BED-format file or, if 'fileName' is not
    given, to STDOUT.

    'compress' can be one of: None, 'gzip', or 'bz2'.
    Compression level is always most compressed/slowest.

    See binfio.featToBed() for explanation of 'nameField'.
    """
    if fileName is None:
        outFile = sys.stdout
    else:
        if compress is None:
            outFile = open (fileName, 'w')
        elif compress == 'gzip':
            outFile = gzip.open (fileName, 'wb', 9)
        elif compress == 'bz2':
            outFile = bz2.BZ2File (fileName, 'w', compresslevel=9)
        else:
            raise exceptions.InvalidArgumentError, \
                  'Invalid compression type: "%s".' % compress

    for feat in coll:
        print >> outFile, feat.toBed (numFields, nameField=nameField)

