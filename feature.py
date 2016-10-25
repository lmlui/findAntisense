#!/bin/env python
"""
Objects representing genomic features and their collections.

...The emphasis is on constructing minimal objects to which you can add
features as you go along, since it is easier for parsing ungodly formats,
recomputing stuff, etc...

General conventions:
  - 'name' is human-readable display name (also for debug out), 'id' is
     unique ID.
  - 1-based coords used for debug output and human readability, 0-based
    half-open (in-between) for internals.

GLOSSARY:
located feature = sequence feature with known coordinates
...stranded versus strandless...
"""


# TODO:
# - RENAME THIS MODULE??  The word 'feature' is too frequently used and module
#   name tends to get rebound -- is there any protection one can use
#   against such a thing?
# - more pydoc!!!!!
# - Exon representation for spliced features should be relative to L coord, not
#   to 5' end??  It might allow me to avoid flipping exons in
#   Spliced.mapToParentRef() and allow me to actually do operations on
#   strandless features (like in Spliced.overlapFeatsToVals()).  [2011/05/05]
# - Should be able to handle negative coords, so that e.g. upstream sequences
#   of a feature could be represented.
# - Maybe have more classes gently adding stuff:
#      LocatedFeature (linear) extends Feature
#      StrandedFeature extends LocatedFeature
#      CircularFeature extends Feature
#   Think of: upstream seqs in both linear and circular nucleic
#   acids, insertion sites, ...
# - Should have a circular feature (e.g. plasmid, circular chromosome) as an
#   extention of the LocatedFeature class.  Must make sure it is a top-level
#   feature (its coords start at 0, it is its own landmark).
# - replace '+'/'-' with constant strings or some sort of method.


import array
import copy
import gzip
import itertools
import operator
import re
import sys

import binfio
import exceptions
import misc


################################################################################
# Classes.
################################################################################


#===============================================================================
# Genomic features.
#===============================================================================


class LocatedFeature (object):
    """A linear feature with known coordinates...

    ...explain why we use 'left' and 'right' instead of 'start' and 'end':
    we see the landmark as being positioned 5'->3', left to right, in a browser
    view, and we describe the endpoints of the sequence relative to the
    landmark, using strand to assign directionality...
    """

    # So many ways to represent strand directionality.
    plusStrandAlphabet, minusStrandAlphabet = set(), set()
    for s in '+PpFfWwTt5': plusStrandAlphabet.add (s)
    for s in '-MmRrCcBb3': minusStrandAlphabet.add (s)
    

    def __init__ (self, landmark, start, end, coordSystem=0, strand=None,
                  strandless=False, computeStrand=False, id=None, name=None,
                  score=None, seq=None, note=None, paranoia=True):
        """...

        TODO: Default coordSystem should be 2, since I use it the most, but
        who knows how much code that will break...  maybe make it a required
        arg for a while, just to see what uses LocatedFeature?

        ...make sure all defaults are clearly mentioned...

        ...define 'landmark'...
        ...'landmark' can be the string name or ID of a landmark, or a reference
        to a feature object...
        TODO: Right now, you are allowed to have a None landmark.
        Callers requiring such an object should ideally use a more generic
        Feature object (i.e. whose locations is not known (e.g. an unplaced
        transcript) or irrelevant (e.g. a chromosome)).

        ...maybe group the description?...

        You cannot use coordSystem={0,1} to represent 0-length features (aka
        junctions), such as insertion sites.  This is due to the existence of
        multiple conventions to represent 0-length features in these coordinate
        systems.  If you want to create a 0-length feature, use 'coordSystem=2',
        which can unambiguously represent these coordinates.

        The value None denotes unknown or unassigned features.

        strandless features are those for which the strand is irrelevant.
        If the strand is important, but not known, set::
            strandless=False, strand=None, computeStrand=False
        
        If no strand is given and the feature is not strandless, the default
        behavior is to compute a strand based on start/end coordinates, i.e. if
        start > end, then we assume minus strand, otherwise plus strand.
        The default can be overridden by setting 'computeStrand=False'.

        Setting 'computeStrand' has no effect if you provide the strand
        explicitly or if the feature is strandless.

        You can specify the coordinate system the feature is in::
            0 = 1-based inclusive (default)
            1 = 0-based inclusive
            2 = 0-based half-open/in-between coords
        (...elaborate above coordinate systems...)
        TODO These should be module constants.
        Internally, we store everything as 0-based in-between coords.

        IDs should be unique to the feature, but names are not.

        'seq'... (TODO Need to check for valid nucs/AAs and convert them to
        a standard form; maybe easier to use BioPython for this?)

        'note' is a string containing any text you want.  It is not used
        for anything, other than serving as a catch-all dump.

        The bulk of this constructor focuses on making sure the arguments
        make sense and that the data is converted into a consistent internal
        representation early on.
        """
        self.paranoia = paranoia  # TODO: this isn't currently used by LocatedFeature, only by Spliced

        # TODO: Check 'landmark' for None and if it is, fall back to some
        # direct reference-to-parent-object mechanism.
        # Or instead, 'landmark' can be a string or another LocatedFeature
        # object, but not None.
        # TODO: strip whitespace from around all these.

        #! Set the easy attributes.
        self.landmark, self.id, self.name, self.seq, self.note \
                       = landmark, id, name, seq, note

        #! Set the score.
        if score is not None:
            try:
                self.score = float (score)
            except ValueError:
                raise exceptions.InvalidArgumentError, \
                      'Value for "score" parameter (%s) is not a number.' \
                      % score
        else:
            self.score = None

        #! Make sure coordinates are numeric.
        try:
            start, end = int (start), int (end)
        except ValueError:
            raise exceptions.InvalidArgumentError, 'Start or end coordinate not an integer.'

        #! Set the strand.
        if strandless:
            if computeStrand:
                raise ArgumentConflictError, \
                      'It is nonsensical to ask to compute a strand for a strandless feature.'
            if strand is not None:
                raise ArgumentConflictError, \
                      'It is nonsensical to provide a strand to a strandless feature.'
        else:
            if computeStrand:
                if strand is not None:
                    raise ArgumentConflictError, \
                          'It is nonsensical to ask to compute a strand when you are already providing a strand.'
                # Figure out the strand based on the coords.
                if start > end: strand = '-'
                else:           strand = '+'
            else:
                if strand is not None:
                    # Standardize strand symbol.
                    if strand in self.plusStrandAlphabet:
                        strand = '+'
                    elif strand in self.minusStrandAlphabet:
                        strand = '-'
                    else:
                        raise exceptions.InvalidArgumentError, 'Invalid strand symbol (%s).' \
                              % strand                    
        self.strand     = strand
        self.strandless = strandless
        
        #! Convert coordinates to 0-based in-between coords for our internal
        #! representation.
        if start > end:
            start, end = end, start
        try:
            coordSystem = int (coordSystem)
        except ValueError:
            raise exceptions.InvalidArgumentError, \
                  'Invalid coordinate system (%s); should be an integer.' \
                  % coordSystem
        if coordSystem == 0:
            # Input coords are 1-based inclusive.
            if start < 1 or end < 1:
                raise exceptions.InvalidArgumentError, \
                      'Invalid coords for coord system 0: [%d,%d].' \
                      % (start, end)
            start -= 1
        elif coordSystem == 1:
            # Input coords are 0-based inclusive.
            if start < 0 or end < 0:
                raise exceptions.InvalidArgumentError, \
                      'Invalid coords for coord system 1: [%d,%d].' \
                      % (start, end)
            end += 1
        elif coordSystem == 2:
            # Input coords are in the same coordinate system as our
            # internal one; just do a safety check.
            if start < 0 or end < 0:
                raise exceptions.InvalidArgumentError, \
                      'Invalid coords for coord system 2: [%d,%d].' \
                      % (start, end)
        else:
            raise exceptions.InvalidArgumentError, \
                  'Invalid coordinate system (%s).' % coordSystem
        self.L, self.R = start, end
        
        if seq is not None:
            if len (seq) != len (self):
               raise exceptions.InvalidArgumentError, \
                   'Sequence length %d of feature "%s" does not equal its coordinate-based length %d.' \
                   % (len (seq), str (self), len (self))


    def __str__ (self):
        """String representation of a LocatedFeature:
        its landmark name, 1-based coordinates, strand (only if stranded), id or name (only if exist).
        """
        if self.id is not None:
            idStr = '{%s}' % self.id
        elif self.name is not None:
            idStr = '{%s}' % self.name
        else:
            idStr = ''
            
        if self.strandless:
            return '%s:%d-%d%s' % (self.landmarkName, self.L+1, self.R, idStr)
        else:
            return '%s:%d-%d(%s)%s' % (self.landmarkName, self.L+1, self.R,
                                     self.strand, idStr)


    def __len__ (self):
        """The length of the feature, in nucleotides."""
        return self.R - self.L


    def getLandmarkName (self):
        """Return landmark ID or name as a string.

        Since the property 'landmark' of a LocatedFeature object can be a
        string name or a reference to another feature object, we need to resolve
        the two ways of getting the landmark name.
        Additionally, we try to find the landmark ID; if it does not exist, we
        try the name, and failing that return the string 'Unnnamed'.
        If there is no landmark, we return the string 'None'.
        """
        if isinstance (self.landmark, str):
            return self.landmark
        elif self.landmark is None:
            return 'None'
        elif self.landmark.id is not None:
            return self.landmark.id
        elif self.landmark.name is not None:
            return self.landmark.name
        else:
            return 'Unnamed'
    landmarkName = property (fget=getLandmarkName)


    def subseq (self, start, end, relativeToOriginal=False):
        """Get the subsequence of a feature as a new LocatedFeature object.

        If 'relativeToOriginal' is True, then the new feature will be relative
        to the original feature (i.e. original is the parent/landmark of the
        new).  The new feature will always be on the plus strand of the
        original.  If the option is False, the new feature will be in the same
        coordinate system as the original and the strand will be preserved.

        'start' and 'end' are in-between coords relative to the original
        feature.

        TODO: This function should be in the more generic feature class
        (without location) when I write it.
        """
        # TODO: should raise custom exceptions, but for now just fubar.
        assert start <= end
        sub = copy.copy (self)
        if sub.seq is not None:
            sub.seq = sub.seq[start:end]
        if relativeToOriginal:
            sub.landmark, sub.L, sub.R, sub.strand = self, start, end, '+'
        else:
            # TODO: UNTESTED!!!
            assert not self.strandless, "TODO: Don't know how to do strandless yet!"
            if self.strand == '+':
                sub.L = self.L + start
                sub.R = self.L + end
            else:
                assert self.strand == '-'
                sub.L = self.R - end
                sub.R = self.R - start
        return sub


    def getAllSubseqs (self, length, relativeToOriginal=False):
        """Get all subsequences of length 'length' and return them as a
        LocatedFeatureCollection, where each subsequence feature has the
        original feature as the parent/landmark.

        TODO: This function should be in the more generic feature class
        (without location) when I write it.
        """

        # TODO: These should raise custom exceptions.  For now, just fubar.
        assert length <= len (self), 'Cannot get a subseq longer than the actual seq.'
        assert self.seq is not None, 'Cannot get a subseq when there is sequence in the parent.'

        subseqs = LocatedFeatureCollection()
        for start in xrange (len (self) - length + 1):
            end = start + length
            subseqs.addFeature (self.subseq (start, end, relativeToOriginal))
        return subseqs


    def mapToParentRef (self, parentCollection=None, disableClassCheck=False):
        """Update the coordinates and landmark/reference of this feature
        so that it is relative to the landmark/reference of its parent
        feature.

        Returns this feature, for convenience.

        Parent of feature A == feature B that is the landmark of feature A.

        Given feature A which uses landmark B, we look up feature B and its
        landmark C.
        Then, we update A so that its landmark is C and its coordinates are
        relative to that landmark.

        If the feature stores a reference to its landmark object in
        'self.landmark', we use that landmark as the parent; otherwise, we look
        up the landmark in 'parentCollection' using the string ID of landmark
        stored in 'self.landmark'.
        So, to use without searching through a parent set, just set feature's
        landmark to parent feature reference.

        (TODO: Get rid of using 'parentCollection',
        force all features to have a direct reference to the landmark object;
        it is too messy this way.)

        Basically, we make it so that the coordinates of feature A are now
        relative to landmark C, instead of landmark B.

        Example of use: let's say that 'self' contains box features of
        snoRNAs.  The coordinates of the box features are relative to their
        respective snoRNA, i.e. their snoRNAs are their landmarks.  We want
        to covert box feature coordinates to be relative to the chromosome.
        We pass in 'parentCollection' containing snoRNAs with coordinates
        relative to the chromosome.  For each box feature, this function will
        find its parent snoRNA in 'parentCollection' and use the snoRNAs
        chromosome coordinates to update the box feature's coordinates so
        that they are on the chromosome.

        If a feature and its parent are both stranded, then the feature strand
        is updated so that it is in the correct orientation relative to the
        parent's landmark.

        If the parent feature is a Spliced object, the algorithm proceeds in
        the same way as if it was a LocatedFeature object (that is, parent exons
        are not considered at all).
        This means that the bounds of this feature could wind up in the introns
        of the parent.
        If this feature is known to be derived only from the exonic regions of
        the parent, do not use this function.
        Instead, make this feature a Spliced feature as well and make it have
        one exon that spans the whole feature.
        The Spliced.mapToParentRef() method will then take care to split the
        exons and compute bounds in such a way that this feature is only over
        the exonic parts of the parent.
        """
        # Get a ref to the parent.
        if isinstance (self.landmark, LocatedFeature):
            # Landmark is already a reference to a feature object.
            parent = self.landmark
        else:
            # Landmark is a string storing the ID of a feature;
            # need to get the object.
            assert parentCollection is not None
            parent = parentCollection.getFeatureByID (self.landmark)
            if parent is None:
                raise exceptions.InvalidArgumentError, \
                      'No parent "%s" of feature: %s' \
                      % (self.landmark, str (self))

        #if not disableClassCheck and isinstance (parent, Spliced):
        #    raise exceptions.InvalidArgumentError, \
        #          'Cannot map non-spliced features (type LocatedFeature) to landmarks of spliced parents (type Spliced), because they may become spliced themselves and then cannot be a LocatedFeature object anymore.  We cannot change the object type.  So, if you want to map a feature to landmarks of spliced parents, make it a Spliced feature, not a LocatedFeature.  If nothing is known about its exons, or the whole feature is exonic, just set it to have one exon spanning the whole feature.'

        # Update the coordinates to the parent's system.
        if self.strandless or parent.strandless:  # TODO
            raise exceptions.FeatureNotImplementedError, \
                  'OOPS! Not sure how to do strandless features yet.'
        self.landmark = parent.landmark
        if parent.strand == '+':
            self.L += parent.L
            self.R += parent.L
            # Leave feature (child) strand the same.
        else:
            assert parent.strand == '-'
            #print 'feature was: %s' % str (self)  #D!!!
            #print 'parent is: %s' % str (parent)  #D!!!
            newL = parent.R - self.R
            newR = parent.R - self.L
            self.L, self.R = newL, newR
            assert self.L <= self.R
            # Feature strand is relative to parent, and the parent's
            # direction is opposite that of the parent's landmark;
            # make direction of feature be relative to parent's landmark
            # by flipping it.
            self.strand = '+' if self.strand == '-' else '-'
            #print 'now feature is: %s' % str (self)  #D!!!
            
        return self
    
    
    def mapToSelf (self):
        """Change coordinates of feature to be relative to self; that is,
        the feature becomes its own landmark.
        
        The resulting interbase coords are [0,L] where L is the feature length.
        
        The resulting strand is always the plus strand.
        
        If the feature has an id, it is used as the landmark; if not, then
        the name is used.  If neither is present, the landmark is None.

        Returns itself, for convenience.
        """
        if self.id is not None:
            self.landmark = self.id
        elif self.name is not None:
            self.landmark = self.name
        else:
            self.landmark = None
        
        self.R      = len (self)
        self.L      = 0
        self.strand = '+'

        return self
    
    
    def toBed (self, numFields, nameField=None):
        """Convert feature to a BED format string with 'numFields' fields
        (columns).

        This method is just a thin wrapper around L{binfio.featToBed},
        so see the docs for that to understand the parameters.
        This also means this method can be used by derived classes like
        L{Spliced}.

        The only reason this method is here is because it seems like good
        practice to bind methods for returning various forms of data
        representation to their objects, rather than force caller to use a
        function in an external module.
        Also this is shorter::
            str = feat.toBed (6)
        than this::
            str = binfio.featToBed (feat, 6)
        """
        return binfio.featToBed (self, numFields, nameField=nameField)


    def overlapFeatsToValsPerPos (self, feats, vals, ignoreStrand=False):
        """Returns a list of values-per-position of 'self' that result from
        covering it with 'feats', each scoring a value from 'vals'.

        Positions for which there is no data get a base-line value of 0.0.

        'vals' are values per feature in 'feats'.
        Although 'feats' and 'vals' can be arbitrary iterables (over
        LocatedFeature objects and numbers, respectively), the iteration must
        return objects in a 1-to-1 mapping, so we can associate each feature
        with its value.

        Features in 'feats' must have the same landmark/reference sequence as
        this feature (at least until overlapFeats() can do a smarter comparison
        between feats, TODO).

        @return: a list of float values, one per position in this feature,
        ordered from L to R coordinates.

        @todo: Use Kevin Karplus' optimization on valsPerPos, like I did in
        the Spliced version of this func.  Also change computeCoverage() and
        other calling funcs.
        """
        valsPerPos = [0.0 for i in xrange (len (self))]
        for feat, val in itertools.izip (feats, vals):
            if overlapFeats (self, feat, ignoreStrand):
                for offset in xrange (max (feat.L, self.L) - self.L,
                                      min (feat.R, self.R) - self.L):
                    valsPerPos[offset] += val
        return valsPerPos


    def computeCoverage (self, otherFeats, ignoreStrand=False):
        """Compute mean per-position coverage of this feature by other
        features.

        'otherFeats' is any iterable yielding feature objects.
        They must have the same landmark/reference sequence as this feature
        (at least until overlapFeats() can do a smarter comparison
        between feats, TODO)
        
        You may pass in features (in 'otherFeats') that overlap this feature
        only partially or not at all.  Only the overlapping regions will
        contribute to coverage.  If 'ignoreStrand' is set to True, features
        on opposite strands will be allowed to overlap.
        
        Each other feature contributes 1x coverage.
        
        TODO: Pass in a "weight vector" for giving feats other counts.
        """
        coveragePerPos = self.overlapFeatsToValsPerPos (
            otherFeats, itertools.repeat (1), ignoreStrand)
        return sum (coveragePerPos) / float (len (self))


    def computeMmi (self, otherFeats, mapCounts, ignoreStrand=False):
        """Compute multiple mapping index (MMI) for this feature imparted by
        other features and their multiple mapping counts.
        
        ...define MMI formally...
        
        You may pass in features (in 'otherFeats') that overlap this feature
        only partially or not at all.  Only the overlapping regions will
        contribute to MMI.  If 'ignoreStrand' is set to True, features
        on opposite strands will be allowed to overlap.

        MMI is always >= 1 because we normalize per-position MMI only to
        positions with 1x coverage or more (TODO: what if coverage is 0
        everywhere?  Code fails with assert on these cases right now.).
        
        Although 'otherFeats' and 'mapCounts' can be arbitrary iterables (over
        LocatedFeature objects and positive integers, respectively), the
        iteration must return objects in a 1-to-1 mapping, so we can associate
        each feature with its multiple mapping count.
        
        Positions with zero coverage are not considered in the MMI calculation.
        
        Each other feature contributes 1x coverage.
        
        TODO: Pass in a "weight vector" for giving feats other counts.
        """
        coveragePerPos = self.overlapFeatsToValsPerPos (
            otherFeats, itertools.repeat (1), ignoreStrand)
        mapCountPerPos = self.overlapFeatsToValsPerPos (
            otherFeats, mapCounts, ignoreStrand)

        assert len (coveragePerPos) == len (mapCountPerPos)

        mmiPerPos = [float (mapCountPerPos[i]) / float (coveragePerPos[i])
                     if coveragePerPos[i] > 0 else 0.0
                     for i in xrange (len (coveragePerPos))]

        # Compute effective length, which is the count of how many
        # positions have nonzero coverage.
        # Normalizing by effective length, rather than actual length,
        # guarantees that, in cases when coverage is 0 for some positions,
        # we always get MMI >= 1.
        effectiveLen = sum ([1 if m > 0 else 0 for m in mmiPerPos])

        # TODO Gentler; warn caller that coverage is zero.
        assert effectiveLen > 0

        return sum (mmiPerPos) / float (effectiveLen)


    def getLociByThreshold (self, valuesPerPos, threshold):
        """Returns a list of loci where each locus spans the run of
        consecutive positions where values from C{valuesPerPos}
        are >= C{threshold}.

        Tip: if C{valuesPerPos} need to be coverage, it can be easily
        obtained with the L{computeCoverage} function.

        Output features all have the same strand, which is the strand of
        C{self}.

        TODO: WiggleWriter should use this to find runs above threshold.
        """
        assert len (valuesPerPos) == len (self)  # TODO: should be an exception
        
        subLoci    = []
        inSubLocus = False
        offset     = 0

        for val in valuesPerPos:
            if not inSubLocus and val >= threshold:
                # Hit beginning of high-coverage run.
                # Start a new sub-locus.
                inSubLocus = True
                subLocusL  = self.L + offset
            elif inSubLocus and val < threshold:
                # Hit end of a high-coverage run.
                # End the sub-locus.
                inSubLocus = False
                subLoci.append (LocatedFeature (
                    self.landmark, subLocusL, self.L + offset,
                    strand=self.strand, coordSystem=2))
            offset += 1
            assert offset <= len (self)
            
        assert offset == len (self)
        
        if inSubLocus:
            # One last sub-locus to save, which must have run to the end
            # of the parent locus.
            subLoci.append (LocatedFeature (
                self.landmark, subLocusL, self.L + offset,
                strand=self.strand, coordSystem=2))
        
        return subLoci

    # Old names pointing to renamed functions (for backwards compatibility).
    mapToParentLandmarks = mapToParentRef


def reverseExons (exons, featLen):
    """A helper function that takes an exon list and reverses it so that it
    is relative to the other end of the feature.

    A new list is made for returning; input exon list is not modified.

    That is, a list of exons with coords relative to the 5' end will be
    converted to a list of exons with coords relative to the 3' end of the
    feature of length 'featLen', and vice versa.

    See 'Spliced' class documentation for the exon list data structure.

    The function is not a method of 'Spliced' primarily because it would be
    bad to have a method that changed a 'Spliced' object's internal exon
    representation to be from the 'wrong end', i.e. inconsistent with the
    spec.
    This func is intended for use with temporary, external exon lists only.

    NB: This function is called by an external module (binfio), so
    careful if refactoring.
    """
    newExons = []
    for i in xrange (len (exons)):
        newExons.append ((featLen - exons[i][1],
                          featLen - exons[i][0]))
    newExons.reverse()
    return newExons


#===============================================================================
# Collections.
#===============================================================================

        
class LocatedFeatureCollection (object):
    """A collection of LocatedFeatures.

    The 'useId' attribue controls how the collection is loaded and indexed,
    namely whether we should index on the ID and coerce names to IDs when
    loading.
    If a feature has a defined ID (i.e. its ID is not None), it must be unique
    to the collection.
    You cannot add a new feature if something with its ID is already in the
    collection.

    TODO: What this class really needs to be is an interface to a full-fledged
    database, maybe as a wrapper around Pygr?  Should have flags/methods for:
      - indexing (turn off/delay for write efficiency)
      - saving/pickling database state to disk

    TODO: When changing a property (e.g. landmark) of an individual feature,
    the binning by that property will break.  Not sure how to really deal with
    this, other than keep a back-pointer in each feature to the collection that
    it is in.  How does Pygr do it?
    """
    def __init__ (self, useId=False):
        self.useId = useId
            
        # Features binned by landmark.
        # Key:   landmark (string).
        # Value: list of features on that landmark.
        self.featuresBinned = {}

        # Features in order that they are added.
        self.featuresOrdered = []

        # Index certain commonly used keys.
        # TODO: This is rapidly approaching what a database does.
        # So, use a real database engine.
        if useId:
            self.featsById = {}  # id --> unique feature
        self.featsByName   = {}  # id --> list of features with that name
        

    def __iter__ (self):
        """Iterate by sorted landmark name, then by order that features
        were added.
        """
        for l in sorted (self.landmarkNames):
            for f in self.featuresBinned[l]:
                yield f


    def __str__ (self):
        return ';  '.join ([str (o) for o in self])


    def __len__ (self):
        return len (self.featuresOrdered)


    def add (self, feature):
        """Add a feature to the collection.

        Returns the added feature, for convenience.

        TODO: Ensure that added feature is of the LocatedFeature type.
        """
        if self.useId and feature.id is not None:
            if self.featsById.has_key (feature.id):
                # TODO: Custom exception.
                assert False, \
                       'Aleady have a feature ID "%s" in this collection.'
            else:
                self.featsById[feature.id] = feature
                
        if feature.name is not None:
            if not self.featsByName.has_key (feature.name):
                self.featsByName[feature.name] = []
            self.featsByName[feature.name].append (feature)

        l = feature.landmarkName
        if not self.featuresBinned.has_key (l):
            self.featuresBinned[l] = []
        self.featuresBinned[l].append (feature)
        
        self.featuresOrdered.append (feature)
        
        return feature

    
    def addFromBed (self, fileName=None, trimWhitespace=True, verbose=0):
        """Add features to this collection, loaded from BED data.

        If 'fileName' is None, we read from STDIN.
        """
        for bedLine in binfio.iterBed (
                fileName, numCols=None, trimWhitespace=trimWhitespace,
                verbose=verbose):
            if self.useId:
                # Try using BED 'name' field for the ID.
                feat = LocatedFeature (
                    bedLine['chrom'], bedLine['chromStart'], bedLine['chromEnd'],
                    coordSystem=2, strand=bedLine['strand'], id=bedLine['name'],
                    score=bedLine['score'])
            else:
                feat = LocatedFeature (
                    bedLine['chrom'], bedLine['chromStart'], bedLine['chromEnd'],
                    coordSystem=2, strand=bedLine['strand'], name=bedLine['name'],
                    score=bedLine['score'])
            self.add (feat)


    def getUniqueNames (self):
        """Returns unique-afied names of features in this collection, in no
        particular order.
        """
        return self.featsByName.keys()
    uniqueNames = property (fget=getUniqueNames)


    def getFeaturesByName (self, name):
        """Get all features with this name as a list, or None if there are no
        such features.
        """
        if self.featsByName.has_key (name):
            return self.featsByName[name]
        else:
            return none


    def getFeatureByID (self, id):
        """Get a feature by its unique identifier, or None if there is no such
        feature.
        """
        if self.featsById.has_key (id):
            return self.featsById[id]
        else:
            return None


    def getLandmarkNames (self):
        """Returns string names of all landmarks used by the features in this
        collection, in no particular order.
        """
        return self.featuresBinned.keys()
    landmarkNames = property (fget=getLandmarkNames)


    def sortByPosition (self, descending=True):
        """Sort the features in this collection by their position, ignoring
        strand info.
        """
        for l in self.landmarkNames:
            self.featuresBinned[l].sort (compare, reverse=descending)
        return self


    def sortBy (self, attributeName, descending=True):
        """Sort the features in this collection by the values of the
        given attribute.
        """
        for l in self.landmarkNames:
            self.featuresBinned[l].sort (
                key=operator.attrgetter (attributeName),
                reverse=descending)
        return self


    def toBedFile (self, fileName=None, compress=None, numFields=6,
                   nameField=None):
        """Write collection to a BED-format file or, if 'fileName' is not
        given, to STDOUT.

        This method is just a thin wrapper around L{binfio.collToBedFile},
        so see the docs for that to understand the parameters.

        The only reason this method is here is because it seems like good
        practice to bind methods for returning various forms of data
        representation to their objects, rather than force caller to use a
        function in an external module.
        Also this is shorter::
            str = coll.toBedFile()
        than this::
            str = binfio.collToBedFile (coll)
        """
        if numFields >= 10:
            raise exceptions.InvalidArgumentError, \
                  "LocatedFeature collection objects do not store splicing info, so can't write BED col 10 or above; you requested numFields = %d." \
                  % (numFields)
        binfio.collToBedFile (
            self, fileName=fileName, compress=compress, numFields=numFields,
            nameField=nameField)
        

    def mapToParentRefs (self, parentCollection):
        """Update the coordinates and landmarks of all features in this
        collection so that they are relative to the landmarks that their
        parent features use.

        Parent of feature A = feature B that is the landmark of feature A.

        Given feature A which uses landmark B, we look up feature B in
        'parentCollection' and its landmark C.  Then, we update A so that
        its landmark is C and its coordinates are in that landmark.

        Basically, we make it so that the coordinates of feature A are now
        relative to landmark C, instead of landmark B, as long as B has C
        as its landmark and B is in 'parentCollection'.

        Example of use: let's say that 'self' contains box features of
        snoRNAs.  The coordinates of the box features are relative to their
        respective snoRNA, i.e. their snoRNAs are their landmarks.  We want
        to covert box feature coordinates to be relative to the chromosome.
        We pass in 'parentCollection' containing snoRNAs with coordinates
        relative to the chromosome.  For each box feature, this function will
        find its parent snoRNA in 'parentCollection' and use the snoRNAs
        chromosome coordinates to update the box feature's coordinates so
        that they are on the chromosome.

        If a feature and its parent are both stranded, then the feature strand
        is updated so that it is in the correct orientation relative to the
        parent's landmark.
        """
        for feature in self:
            feature.mapToParentRef (parentCollection)
            
        # Update the binning, since the landmarks have changed.
        newBins = {}
        for feature in self:
            l = feature.landmarkName
            if not newBins.has_key (l):
                newBins[l] = []
            newBins[l].append (feature)
        self.featuresBinned = newBins
        

    def toGroups (self, spacing=0, stranded=True):
        """Merge overlapping/adjacent features in this collection into groups;
        features can be from different landmarks/strands.

        Features 'spacing' nucleotides apart or closer are merged.  If
        'spacing' is 0, only immediately adjacent features are merged into
        groups.

        'stranded' controls whether we merge overlapping features on
        opposite strands.

        Returns a LocatedFeatreCollection of groups.  A group is a strandless
        LocatedFeature object with a new attribute called 'reads' that is an
        ordered list of reads that went into making that group.

        TODO: rename 'reads' to something less domain-dependent
        """
        ret = LocatedFeatureCollection()
        
        if stranded:
            for featsSameLandmark in self.featuresBinned.values():
                # Make iterator over features on the same strand;
                # use lazy iteration for memory efficiency.
                def featsSameStrand (feats, strand):
                    for f in feats:
                        assert not f.strandless  # TODO gentler
                        if f.strand == strand:
                            yield f

                # Separate features by strand and do each strand separately,
                # because groupsFromContig() does not look at strand when
                # making groups.
                for strand in ('+', '-'):
                    groups, featsByGroup = groupsFromContig (
                        featsSameStrand (featsSameLandmark, strand),
                        spacing=spacing)
                    assert len (groups) == len (featsByGroup)
                    for i in xrange (len (groups)):
                        group = groups[i]
                        # Add list of reads (ordered by coord)
                        # to the group.
                        group.reads = featsByGroup[i]
                        # Add strand info to group (groupsFromContig() does
                        # not do strand info).
                        group.strandless = False
                        group.strand     = strand
                        # Save the group in our new collection.
                        ret.addFeature (group)
                        
        else:
            # Ignore all strand information.
            for featsSameLandmark in self.featuresBinned.values():
                groups, featsByGroup = groupsFromContig (
                    featsSameLandmark, spacing=spacing)
                assert len (groups) == len (featsByGroup)
                for i in xrange (len (groups)):
                    # Add list of reads (ordered by coord)
                    # to the group.
                    groups[i].reads = featsByGroup[i]
                    # Save the group in our new collection.
                    ret.addFeature (groups[i])

        return ret


################################################################################
# Factory functions.
################################################################################
          
def bedToLocusColl (fileName=None, useId=False, verbose=0, trimWhitespace=True):
    """Factory function for creating a collection from features in
    BED format.

    If no filename is given, we read from STDIN.

    TODO: Shouldn't store landmark as a string.  Should automatically create
    a 'coordinate-less' top-level feature object to represent, to get away from
    this craziness of storing stuff as strings.
    """
    retColl = LocatedFeatureCollection (useId=useId)
    retColl.addFromBed (
        fileName, trimWhitespace=trimWhitespace, verbose=verbose)
    return retColl
    
# Alias for backward compatiblity.
locatedFeatureCollectionFromBED = bedToLocusColl


################################################################################
# Other module functions.
################################################################################

def groupsFromContig (features, spacing=0):
    """Merge overlapping/adjacent features into numbered groups, with the assumption that
    all features are on same contig and strand (or that strand doesn't matter).

    Group numbers are integers beginning at 0.

    @attention: The strand of input features is never checked by this function!  It is up to the CALLER whether
        to filter features by strand before passing or to allow this function to group together features on opposite strands.

    @param spacing: Features C{spacing} bases apart or closer are merged into a group.
    @type spacing: int

    @param features: Any iterable that yields L{LocatedFeature} objects.  Input features are not modified.  They do not have to be sorted.
    @type features: iterable

    @rtype: The tuple C{(groups, featsByGroup)}
    @return: C{groups} is a list of groups (strandless L{LocatedFeature} objects) indexed by
        the integer group number.
        Returned groups are strandless; the caller must assign strand, if that is
        desired.  However, they have a landmark assigned (taken from input features).
        C{featsByGroup} is a dict where the key is the group number and the value is a list
        of features in that group, ordered by their location.
    """
    # TODO: try/except to check for invalid params
    spacing = int (spacing)
    assert spacing >= 0  # TODO: gentler, raise exception
    
    # Ordered by location.  Indexed by groupNum.
    # Each group is a strandless LocatedFeature object.
    groups = []
    # Key: groupNum; value: list of features in group.
    featsByGroup = dict()

    # Gets incremented so that offset begins at 0.
    groupNum = -1
    # Init group coords out of bound (i.e. no overlap with any feature
    # possible) to denote that no groups have been started yet.
    groupL = -1 - spacing
    groupR = groupL

    landmark = None  # for safety checks

    # Identify groups by going through the sorted feature list and
    # merging features into a current group if they overlap.
    for feat in sorted (features, compare):
        if landmark is not None:
            # Feats must all be on same landmark.
            assert feat.landmark == landmark
            # TODO: do above gentler (raise BadInputData exception or some such??)
        else:
            # Must be on first feature.
            landmark = feat.landmark

        if feat.L <= (groupR + spacing):
            # Feature overlaps with or is sufficiently adjacent to current
            # group.  "Merge" the feature into the group.
            groupR = max (feat.R, groupR)
            assert feat.L >= groupL
        else:
            # Feature is disjoint from group.
            if groupL >= 0:
                # Previous group exists (groupL < 0 iff no groups yet started);
                # save it.
                assert groupNum >= 0
                groups.append (LocatedFeature (
                    landmark, groupL, groupR, strandless=True, coordSystem=2))
                
            # Start a new group using the feature's bounds.
            groupNum += 1
            groupL, groupR = feat.L, feat.R

        # Save feature to the current group.
        assert groupNum >= 0
        if not featsByGroup.has_key (groupNum):
            featsByGroup[groupNum] = []
        featsByGroup[groupNum].append (feat)

    if groupL >= 0:
        # One last group left on this strand; save it.
        groups.append (LocatedFeature (
            landmark, groupL, groupR, strandless=True, coordSystem=2))
    else:
        # Special case; occurs when we have zero features
        # as input, so no groups made.
        assert groupNum           == -1
        assert len (groups)       ==  0
        assert len (featsByGroup) ==  0

    return (groups, featsByGroup)

